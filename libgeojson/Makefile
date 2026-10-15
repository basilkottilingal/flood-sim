DEPS :=
include ../lib.mk

lxr: json.l
	lexer -d -o lexer.h json.l  ## -d is for debug
parser: json.y
	bison -d json.y 

# [uncomment this] In case of rebuilding the parser 
# ast.o : lxr parser


