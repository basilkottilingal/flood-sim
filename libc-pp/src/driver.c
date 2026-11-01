/* cwrap: a gcc-compatible front end. For each .c file:
 *
 *     gcc -E  ->  translate_file()  ->  gcc -x cpp-output
 *
 * then links as gcc would. Linemarkers survive the middle step, so diagnostics
 * and debug info point at the original source.
 *
 *   CWRAP_CC       compiler to run (default gcc)
 *   CWRAP_VERBOSE  print each command
 *   CWRAP_KEEP     keep the intermediate directory
 */
#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

#include "translate.h"

enum stage {
	ST_ALL,		/* every gcc run */
	ST_PP,		/* preprocessing only */
	ST_LINK,	/* link only; also linker inputs, in command-line order */
};

struct item {
	enum stage stage;
	const char *flag, *value;	/* value is a separate argument */
	const char *src;		/* .c input */
	char *obj;			/* its object, when linking */
};

struct vec {
	const char **v;
	size_t n, cap;
};

static const struct { const char *name; enum stage stage; } value_flags[] = {
	{ "-I", ST_PP }, { "-D", ST_PP }, { "-U", ST_PP }, { "-A", ST_PP },
	{ "-include", ST_PP }, { "-imacros", ST_PP }, { "-isystem", ST_PP },
	{ "-iquote", ST_PP }, { "-idirafter", ST_PP }, { "-iprefix", ST_PP },
	{ "-iwithprefix", ST_PP }, { "-Xpreprocessor", ST_PP },
	{ "-MF", ST_PP }, { "-MT", ST_PP }, { "-MQ", ST_PP },
	{ "-L", ST_LINK }, { "-l", ST_LINK }, { "-T", ST_LINK },
	{ "-u", ST_LINK }, { "-z", ST_LINK }, { "-Xlinker", ST_LINK },
	{ "-Xassembler", ST_ALL }, { "-B", ST_ALL },
};

/* These may also be joined to their value, as in -Ifoo. */
static const char joined_flags[][3] = {
	"-I", "-D", "-U", "-A", "-L", "-l", "-T", "-z", "-B",
};

static const char *const pp_flags[] = {
	"-MD", "-MMD", "-MP", "-nostdinc", "-undef",
};

/* Modes that need no translation: hand the whole command line to gcc. */
static const char *const passthrough_flags[] = {
	"-E", "-M", "-MM", "-x",
};

static const char *cc;
static char tmpdir[4096];
static int verbose;

static void die(const char *fmt, ...)
{
	va_list ap;

	fputs("cwrap: ", stderr);
	va_start(ap, fmt);
	vfprintf(stderr, fmt, ap);
	va_end(ap);
	fputc('\n', stderr);
	exit(1);
}

static void push(struct vec *v, const char *s)
{
	if (v->n + 2 > v->cap) {
		v->cap = v->cap ? v->cap * 2 : 32;
		v->v = realloc(v->v, v->cap * sizeof(*v->v));
		if (!v->v)
			die("out of memory");
	}
	v->v[v->n++] = s;
	v->v[v->n] = NULL;
}

static char *fmt(const char *format, ...)
{
	va_list ap;
	char *s;

	va_start(ap, format);
	if (vasprintf(&s, format, ap) < 0)
		die("out of memory");
	va_end(ap);
	return s;
}

static int in_list(const char *s, const char *const *list, size_t n)
{
	for (size_t i = 0; i < n; i++)
		if (strcmp(s, list[i]) == 0)
			return 1;
	return 0;
}

#define IN_LIST(s, list) in_list(s, list, sizeof(list) / sizeof(*(list)))

static const char *extension(const char *path)
{
	const char *base = strrchr(path, '/');
	const char *dot = strrchr(base ? base + 1 : path, '.');

	return dot ? dot + 1 : "";
}

/* "dir/foo.c" -> "foo" */
static char *stem(const char *path)
{
	const char *base = strrchr(path, '/');
	char *s = strdup(base ? base + 1 : path);
	char *dot = strrchr(s, '.');

	if (dot)
		*dot = '\0';
	return s;
}

static int run(struct vec *cmd)
{
	pid_t pid;
	int status;

	if (verbose) {
		for (size_t i = 0; i < cmd->n; i++)
			fprintf(stderr, "%s%s", i ? " " : "+ ", cmd->v[i]);
		fputc('\n', stderr);
	}
	pid = fork();
	if (pid < 0)
		die("fork: %s", strerror(errno));
	if (pid == 0) {
		execvp(cmd->v[0], (char *const *)cmd->v);
		fprintf(stderr, "cwrap: %s: %s\n", cmd->v[0], strerror(errno));
		_exit(127);
	}
	if (waitpid(pid, &status, 0) < 0)
		die("waitpid: %s", strerror(errno));
	if (WIFEXITED(status))
		return WEXITSTATUS(status);
	return 128 + (WIFSIGNALED(status) ? WTERMSIG(status) : 0);
}

static void cleanup(void)
{
	struct vec rm = { 0 };

	if (!tmpdir[0])
		return;
	if (getenv("CWRAP_KEEP")) {
		fprintf(stderr, "cwrap: kept %s\n", tmpdir);
		return;
	}
	verbose = 0;
	push(&rm, "rm");
	push(&rm, "-rf");
	push(&rm, tmpdir);
	run(&rm);
}

static void add_stage_flags(struct vec *cmd, const struct item *items,
			    size_t n, enum stage stage)
{
	for (size_t i = 0; i < n; i++) {
		if (items[i].src || items[i].stage != stage)
			continue;
		push(cmd, items[i].flag);
		if (items[i].value)
			push(cmd, items[i].value);
	}
}

static void make_tmpdir(void)
{
	const char *base = getenv("TMPDIR");

	snprintf(tmpdir, sizeof(tmpdir), "%s/cwrap.XXXXXX",
		 base && *base ? base : "/tmp");
	if (!mkdtemp(tmpdir)) {
		int err = errno;

		tmpdir[0] = '\0';
		die("mkdtemp: %s", strerror(err));
	}
}

int main(int argc, char **argv)
{
	struct item *items = calloc((size_t)argc, sizeof(*items));
	size_t nitems = 0, nsrc = 0;
	const char *output = NULL, *dep_file = NULL;
	int compile_only = 0, assemble_only = 0, passthrough = 0;
	int want_deps = 0, have_dep_target = 0;
	int other_sources = 0, nobj = 0;

	if (!items)
		die("out of memory");
	cc = getenv("CWRAP_CC");
	if (!cc || !*cc)
		cc = "gcc";
	verbose = getenv("CWRAP_VERBOSE") != NULL;

	for (int i = 1; i < argc; i++) {
		const char *a = argv[i];
		struct item *it = &items[nitems];
		size_t j;

		if (a[0] != '-' || a[1] == '\0') {
			const char *ext = extension(a);
			static const char *const sources[] = {
				"s", "S", "i", "cc", "cpp", "cxx", "C", "m",
			};

			if (strcmp(a, "-") == 0) {
				passthrough = 1;
			} else if (strcmp(ext, "c") == 0) {
				it->src = a;
				nsrc++;
				nitems++;
				continue;
			}
			other_sources |= IN_LIST(ext, sources);
			it->stage = ST_LINK;
			it->flag = a;
			nitems++;
			continue;
		}
		if (strcmp(a, "-c") == 0) {
			compile_only = 1;
			continue;
		}
		if (strcmp(a, "-S") == 0) {
			assemble_only = 1;
			continue;
		}
		if (strcmp(a, "-o") == 0) {
			if (++i == argc)
				die("missing filename after -o");
			output = argv[i];
			continue;
		}
		if (IN_LIST(a, passthrough_flags))
			passthrough = 1;
		if (strcmp(a, "-MD") == 0 || strcmp(a, "-MMD") == 0)
			want_deps = 1;

		it->flag = a;
		it->stage = ST_ALL;
		nitems++;

		if (IN_LIST(a, pp_flags)) {
			it->stage = ST_PP;
			continue;
		}
		if (strncmp(a, "-Wp,", 4) == 0) {
			it->stage = ST_PP;
			continue;
		}
		if (strncmp(a, "-Wl,", 4) == 0) {
			it->stage = ST_LINK;
			continue;
		}
		for (j = 0; j < sizeof(value_flags) / sizeof(*value_flags); j++) {
			if (strcmp(a, value_flags[j].name) != 0)
				continue;
			if (++i == argc)
				die("missing argument to %s", a);
			it->stage = value_flags[j].stage;
			it->value = argv[i];
			if (strcmp(a, "-MF") == 0)
				dep_file = argv[i];
			if (strcmp(a, "-MT") == 0 || strcmp(a, "-MQ") == 0)
				have_dep_target = 1;
			break;
		}
		if (it->value)
			continue;
		for (j = 0; j < sizeof(joined_flags) / sizeof(*joined_flags); j++) {
			size_t len = strlen(joined_flags[j]);
			size_t k;

			if (strncmp(a, joined_flags[j], len) != 0)
				continue;
			for (k = 0; k < sizeof(value_flags) / sizeof(*value_flags); k++)
				if (strcmp(value_flags[k].name, joined_flags[j]) == 0)
					it->stage = value_flags[k].stage;
			break;
		}
		if (strncmp(a, "-MF", 3) == 0) {
			it->stage = ST_PP;
			dep_file = a + 3;
		} else if (strncmp(a, "-MT", 3) == 0 || strncmp(a, "-MQ", 3) == 0) {
			it->stage = ST_PP;
			have_dep_target = 1;
		}
	}

	if (passthrough || nsrc == 0) {
		argv[0] = (char *)cc;
		execvp(cc, argv);
		die("%s: %s", cc, strerror(errno));
	}
	if ((compile_only || assemble_only) && output && nsrc > 1)
		die("-o cannot be used with -c or -S and multiple files");
	if ((compile_only || assemble_only) && other_sources)
		die("mixing .c with other sources under -c or -S is not supported");

	make_tmpdir();
	atexit(cleanup);

	for (size_t i = 0; i < nitems; i++) {
		struct item *it = &items[i];
		struct vec pp = { 0 }, comp = { 0 };
		const char *result;
		char *s, *pre, *post;
		int rc;

		if (!it->src)
			continue;
		s = stem(it->src);
		pre = fmt("%s/%d-%s.i", tmpdir, nobj, s);
		post = fmt("%s/%d-%s.t.i", tmpdir, nobj, s);
		if (compile_only || assemble_only)
			result = output ? output
					: fmt("%s.%s", s, assemble_only ? "s" : "o");
		else
			result = it->obj = fmt("%s/%d-%s.o", tmpdir, nobj, s);
		nobj++;

		push(&pp, cc);
		push(&pp, "-E");
		add_stage_flags(&pp, items, nitems, ST_ALL);
		add_stage_flags(&pp, items, nitems, ST_PP);
		if (want_deps && !dep_file) {
			char *base = strdup(result);
			char *dot = strrchr(base, '.');

			if (dot)
				*dot = '\0';
			push(&pp, "-MF");
			push(&pp, fmt("%s.d", base));
		}
		if (want_deps && !have_dep_target) {
			push(&pp, "-MT");
			push(&pp, compile_only || assemble_only
				  ? result : fmt("%s.o", s));
		}
		push(&pp, it->src);
		push(&pp, "-o");
		push(&pp, pre);
		rc = run(&pp);
		if (rc)
			return rc;

		if (translate_file(pre, post))
			return 1;

		push(&comp, cc);
		add_stage_flags(&comp, items, nitems, ST_ALL);
		push(&comp, assemble_only ? "-S" : "-c");
		push(&comp, "-x");
		push(&comp, "cpp-output");
		push(&comp, post);
		push(&comp, "-o");
		push(&comp, result);
		rc = run(&comp);
		if (rc)
			return rc;
	}

	if (compile_only || assemble_only)
		return 0;

	{
		struct vec link = { 0 };

		push(&link, cc);
		add_stage_flags(&link, items, nitems, ST_ALL);
		for (size_t i = 0; i < nitems; i++) {
			if (items[i].src) {
				push(&link, items[i].obj);
			} else if (items[i].stage == ST_LINK) {
				push(&link, items[i].flag);
				if (items[i].value)
					push(&link, items[i].value);
			}
		}
		if (output) {
			push(&link, "-o");
			push(&link, output);
		}
		return run(&link);
	}
}
