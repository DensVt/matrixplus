# Определение операционной системы
OS := $(shell uname)

# Компилятор и флаги
CC = g++
SRC = s21_matrix_oop.cc
OBJ_OOP_O = s21_matrix_oop.o
OBJ_TEST_O = s21_matrix_test.o
ARCHIVER_OOP = s21_matrix_oop.a
ARCHIVER_TEST = s21_matrix_test.a
FLAGS= -Wall -Werror -Wextra -std=c++17
TESTFLAGS=-lgtest

# Выбор команды для открытия HTML-отчета
ifeq ($(OS),Linux)
    OPEN_CMD = xdg-open report/index.html
endif
ifeq ($(OS),Darwin)
    OPEN_CMD = open report/index.html
else
    OPEN_CMD = powershell.exe Start-Process 'report\\index.html'
endif

# Цели
all: gcov_report

gcov_report: clean s21_matrix_test.a
	$(CC) -coverage $(SRC) $(ARCHIVER_TEST) -o s21_matrix_oop.out $(TESTFLAGS) -lpthread
	./s21_matrix_oop.out
	lcov -t "s21_matrix_oop" -o s21_matrix_oop.info -c -d ./ --no-external
	genhtml -o report s21_matrix_oop.info
	$(OPEN_CMD)

test: s21_matrix_oop.a
	$(CC) $(FLAGS) s21_matrix_test.cc $(ARCHIVER_OOP) -o s21_matrix_oop.out $(TESTFLAGS)
	./s21_matrix_oop.out

s21_matrix_oop.a:
	$(CC) -c $(FLAGS) $(SRC) -o $(OBJ_OOP_O)
	ar -rcs $(ARCHIVER_OOP) $(OBJ_OOP_O)

s21_matrix_test.a:
	$(CC) -c $(FLAGS) s21_matrix_test.cc -o $(OBJ_TEST_O)
	ar -rcs $(ARCHIVER_TEST) $(OBJ_TEST_O)

clean:
	rm -rf *.out *.o *.gcno *.info s21_matrix_oop.a *.gcda *.a *.o
	-rm -rf report

# Условия для команды leaks
ifeq ($(OS), Darwin)
leaks:
	leaks --atExit -- ./s21_matrix_oop.out
endif

ifeq ($(OS), Linux)
leaks:
	valgrind --leak-check=full ./s21_matrix_oop.out
endif

# Цель для clang-format
clang:
	cp ../materials/linters/.clang-format .
	clang-format -n *.cc
	clang-format -n *.h
	clang-format -i *.cc
	clang-format -i *.h
	rm .clang-format

.PHONY: clean clang all leaks