#ifndef _GEOJSON_PARSER_AST_H_
#define _GEOJSON_PARSER_AST_H_

  #include <stdlib.h>

  typedef struct Ast
  {
    struct Ast * child;
    struct Ast * parent;
    struct Ast * next;
    int type;
    int yylineno;
  } Ast;

  typedef struct
  {
    Ast a;
    char * token;
  } AstLeaf;

  typedef struct Allocator
  {
    char * mem;
    struct Allocator * next;
    size_t left;
  } Allocator;

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
  void * allocate (size_t s);
  void   ast_free_all ();
  Ast *  ast_json   (Ast * val);
  Ast *  ast_val    (Ast * );
  Ast *  ast_obj    (Ast * members);
  Ast *  ast_pair   (Ast * str, Ast * val);
  Ast *  ast_array  (Ast * elements);
  Ast *  ast_append (Ast *, Ast *);
  Ast *  ast_leaf   (char * yyText,int yyLen, int tkn);
#endif
