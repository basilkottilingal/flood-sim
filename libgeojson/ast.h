#ifndef _GEOJSON_PARSER_AST_H_
#define _GEOJSON_PARSER_AST_H_

  typedef struct Ast
  {
    struct Ast ** children, * parent;
  } Ast;

  typedef struct
  {
    Ast ** children, * parent;
    char * start, * end;
  } AstLeaf;

  /* APIs */
  void  ast_free_all ();
  Ast * ast1         (Ast *);
  Ast * ast2         (Ast *, Ast *);
  Ast * ast3         (Ast *, Ast *, Ast *);
  Ast * ast_leaf     (char *);
#endif
