CC = g++
CFLAGS = --std=c++17 -Wall -Wextra -Werror
LIBS = -lgtest -lgtest_main -lpthread

LIBRARY = s21_matrix_oop.a

TEST_SRC = $(wildcard tst_src/*.cpp)
TEST_OBJ = $(TEST_SRC:.cpp=.o)

LIB_SRC = $(wildcard lib_src/*.cpp)
LIB_OBJ = $(LIB_SRC:.cpp=.o)

HEADERS = $(wildcard *.h) $(wildcard lib_src/*.h) $(wildcard tst_src*.h)

PROGRAM = test

all: test clean

test: $(LIBRARY) $(TEST_OBJ)
	@$(CC) $(CFLAGS) $(TEST_OBJ) $(LIBRARY) -o $(PROGRAM) $(LIBS)
	@./$(PROGRAM) || true

gcov_report: $(TEST_OBJ)
	@mkdir report
	@$(CC) $(CFLAGS) $(TEST_OBJ) $(LIB_SRC) --coverage -o tests $(LIBS)
	@./tests
	@lcov -t "tests" -o test.info -c -d .
	@genhtml -o report test.info
	@open report/index.html

s21_matrix_oop.a: $(LIB_OBJ)
	@ar rcs $(LIBRARY) $(LIB_OBJ)
	@ranlib $(LIBRARY)

%.o: %.cpp
	@$(CC) $(CFLAGS) -c $< -o $@

check: $(TEST_OBJ) $(LIBRARY)
	@cppcheck $(TEST_SRC) $(LIB_SRC) $(HEADERS)

clang-check:
	@cp ../materials/linters/.clang-format .
	@clang-format -n -style=Google $(TEST_SRC) $(LIB_SRC) $(HEADERS)

	@rm -f .clang-format
	
clang-format:
	@cp ../materials/linters/.clang-format .
	@clang-format -i -style=Google $(TEST_SRC) $(LIB_SRC) $(HEADERS)

	@rm -f .clang-format

valgrind: test
	valgrind --tool=memcheck --leak-check=yes --leak-check=full --show-leak-kinds=all ./test > /dev/null

clean:
	@rm -rf $(TEST_OBJ) $(LIB_OBJ) $(LIBRARY) $(PROGRAM) *.css *.html *.gcno *.gcda *.info tests report *.dSYM