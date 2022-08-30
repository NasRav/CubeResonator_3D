NAME		=	CubeResonator_3D

CLASSES		=	Resonator Explicit
FUNCTIONS	=	CubeResonator_3D $(CLASSES)
SRCS		=	$(patsubst %,%.cpp,$(FUNCTIONS))
OBJS		=	$(patsubst %,%.o,$(FUNCTIONS))
HEADERS		=	$(patsubst %,%.h,$(CLASSES))

CPP			=	g++
FLAGS		=	-Wall -Wextra -Werror -fopenmp
RM			=	/bin/rm -f

.PHONY: all clean fclean re

all:		$(NAME)

$(NAME):	$(OBJS) $(HEADERS)
			$(CPP) $(OBJS) $(FLAGS) -o $(NAME)

clean:
			$(RM) $(OBJS)

fclean:		clean
			$(RM) $(NAME)

re:			fclean	all
