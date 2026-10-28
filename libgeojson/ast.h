#ifndef _GEOJSON_PARSER_AST_H_
#define _GEOJSON_PARSER_AST_H_

  #include <stdint.h>
  #include "liballoc/allocator.h"

  typedef struct Ast
  {
    struct Ast * child;
    struct Ast * parent;
    struct Ast * next;
    int type;
    int n;
  } Ast;

  typedef struct
  {
    Ast a;
    char * token;
    int yylineno;
    int yycolno;
  } AstLeaf;

  typedef struct AstRoot
  {
    Ast a;
    Allocator * alloc;
    uint64_t magic;
  } AstRoot;

  /* in order to avoid circular dependency, we set ast.h as the root
  .. [and thus avoid dependency on lexer/parser definitions]. so we
  .. need to define enum for lexer tokens and parser states.
  */
  #define ast_token(tkn)  ast_kind_##tkn
  typedef enum
  {
    ast_token (jnull) = 1, 
    ast_token (jstring), 
    ast_token (jtrue), 
    ast_token (jfalse), 
    ast_token (jnumber),

    ast_token (jjson) = 256,
    ast_token (jobject),
    ast_token (jvalue),
    ast_token (jpair),
    ast_token (jarray),
  } AstToken;
  #define ast_is_leaf(tkn) (tkn < 256)

  /* APIs */
  Ast *  ast_json   (Allocator *, Ast * val);
  Ast *  ast_val    (Allocator *, Ast * );
  Ast *  ast_obj    (Allocator *, Ast * members);
  Ast *  ast_pair   (Allocator *, Ast * str, Ast * val);
  Ast *  ast_array  (Allocator *, Ast * elements);
  Ast *  ast_append (Allocator *, Ast *, Ast *);
  Ast *  ast_leaf   (Allocator *, char * yyText,int yyLen, int tkn);
  Allocator * ast_allocator (Ast * root);
  void   ast_free   (Ast * root);
  void   ast_print  (Ast * root);
#endif
