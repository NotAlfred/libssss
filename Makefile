##
#  sssslib  -  Copyright 2019 NotAlfred
##

NAME		=       libssss

CC		=	gcc -g

INC		=	includes

CFLAGS		=	-W -Wall -Wextra -pedantic -O2 -std=c99 -I$(INC)#-fPIC

LDFLAGS		=	-lgmp #-shared

RM		=	rm -f

SOURCES		=	src

FILES		=	$(SOURCES)/ssss.c

OBJS		=	$(FILES:.c=.o)

all		:
			@$(MAKE) --no-print-directory $(NAME)

libssss		:	$(OBJS)
			$(CC) $(CFLAGS) $(OBJS) -o $(NAME) $(LDFLAGS)

clean		:
			$(RM) $(OBJS) $(OBJC)

fclean		:	clean
			$(RM) $(NAME)

re		:	fclean all

.PHONY		:	all clean fclean re
