#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>

// Constants
#define PI 3.14159265358979323846
#define E 2.71828182845904523536

typedef enum {
    NODE_NUMBER,
    NODE_VARIABLE,
    NODE_OPERATOR,
    NODE_FUNCTION
} NodeType;

typedef struct Node {
    NodeType type;
    union {
        double number;
        char variable;
        struct {
            char op;
            struct Node *left, *right;
        } op;
        struct {
            char func[10];
            struct Node *child;
            int base; // Base for logarithmic functions (0 if not used)
            char var_base; // If base is a variable
            struct Node *base_node; // If base is a function
        } func;
    };
} Node;

typedef struct {
    double** data;
    int rows;
    int cols;
} Matrix;

// Function declarations
Node* parseExpression(const char **str);
Node* parseTerm(const char **str);
Node* parsePower(const char **str);
Node* parseFunction(const char **str);
Node* parseFactor(const char **str);
void printTree(Node *node);
void freeTree(Node *node);
double evaluate(Node *node, double x);
void printMenu(int* choice);
Node* getFunction(char* input, size_t size);
int compareSigns(double a, double b);
double bisection(Node* root, double start, double end, double epsilon, int maxIterations);
double applyBisection(Node* root);
double regulaFalsi(Node* root, double start, double end, double epsilon, int maxIterations);
double applyRegulaFalsi(Node* root);
double forwardDifferenceDerivative(Node* root, double x, double h);
double backwardDifferenceDerivative(Node* root, double x, double h);
double centralDifferenceDerivative(Node* root, double x, double h);
void applyDifferentiation(Node* root);
double newtonRaphson(Node* root, double x0, double epsilon, int maxIterations);
double applyNewtonRaphson(Node* root);
double simpson1_3(Node* root, double start, double end, double n);
double applySimpson1_3(Node* root);
double simpson3_8(Node* root, double start, double end, double n);
double applySimpson3_8(Node* root);
double* applySimpson(Node* root);
double trapezoidal(Node* root, double start, double end, double n);
double applyTrapezoidal(Node* root);
Matrix createMatrix(int rows, int cols);
void freeMatrix(Matrix m);
Matrix getMatrix();
Matrix getAugmentedMatrix();
Matrix getSquareMatrix();
Matrix multiplyMatrices(Matrix a, Matrix b);
Matrix* separateAugmentedMatrix(Matrix a);
int isSquareMatrix(Matrix a);
int isAugmentedMatrix(Matrix a);
void printAugmentedMatrixEquation(Matrix a, char* title, char var, int indexed);
void printMatrixEquation(Matrix coefs, Matrix results, char* title, char var, int indexed);
Matrix upperTriangularGaussianElimination(Matrix a, double* det, int* rank);
Matrix lowerTriangularGaussianElimination(Matrix a, double* det, int* rank);
double* solveAugmentedMatrix(Matrix a);
double* solveEquationSystem(Matrix coefs, Matrix results);
int isSymmetric(Matrix a);
int isPositiveDefinite(Matrix a);
double sumUpTo(int n);
double factorial(int n);
Matrix getTranspose(Matrix a);
Matrix* choleskyFactorization(Matrix a);
void printMatrixDecomposition(char* title, Matrix A, char* a, Matrix L, char* l, Matrix U, char* u);
double* choleskyDecomposition(Matrix a);
void applyCholesky(Matrix a);
double determinant(Matrix a);
int isInvertible(Matrix a);
Matrix inverseMatrix(Matrix a);
void printMatrix(Matrix m);
Matrix multiplyMatrixWithTranspose(Matrix a);
void swapRows(Matrix* a, int row1, int row2);
double* gaussSeidal(Matrix a, double epsilon, int maxIterations);
double* applyGaussSeidal(Matrix matrix);
double** getForwardDifferenceTable(double* x, double* y, int n);
void printForwardDifferenceTable(double* x, double** arr, int n);
double product(double* x, int size, double searchX, int iteration);
void printProducts(double* x, int size, double searchX, int iteration);
double gregoryNewtonForward(double* x, double* y, int n, double searchX);
double applyGregoryNewtonForward();

// MAIN
int main() {
    char input[256];
    int run = 1;
    int choice, i;

    double rootValue;
    Node* root;

    double det;
    double* results;
    int rank;
    Matrix matrix, temp;
    Matrix* matrices;

    while (run) {
        printMenu(&choice);

        switch (choice) {
            case 0:
                run = 0;
                break;
            case 1:
                root = getFunction(input, sizeof(input));
                // Print function tree
                printf("Function: ");
                printTree(root);
                printf("\n");
                // BISECTION
                rootValue = applyBisection(root);
                // Free memory
                freeTree(root);
                break;
            case 2:
                root = getFunction(input, sizeof(input));
                // Print function tree
                printf("Function: ");
                printTree(root);
                printf("\n");
                // REGULA-FALSI
                rootValue = applyRegulaFalsi(root);
                // Free memory
                freeTree(root);
                break;
            case 3:
                root = getFunction(input, sizeof(input));
                // Print function tree
                printf("Function: ");
                printTree(root);
                printf("\n");
                // NEWTON-RAPHSON
                rootValue = applyNewtonRaphson(root);
                // Free memory
                freeTree(root);
                break;
            case 4:
                matrix = getSquareMatrix();
                // Print Matrix
                printf("\nMatrix:\n");
                printMatrix(matrix);
                // INVERSE OF SQUARE MATRIX
                temp = inverseMatrix(matrix);
                // Print Inverse Matrix
                if (temp.data != NULL) {
                    printf("Inverse Matrix:\n");
                    printMatrix(temp);
                }
                else {
                    printf("\n");
                }
                // Free memory
                freeMatrix(matrix);
                freeMatrix(temp);
                break;
            case 5:
                // Get augmented matrix
                matrix = getAugmentedMatrix();
                // Print Matrix
                printf("\nMatrix:\n");
                printMatrix(matrix);
                // CHOLESKY DECOMPOSITION
                applyCholesky(matrix);
                // Free memory
                freeMatrix(matrix);
                break;
            case 6:
                // Get augmented matrix
                matrix = getAugmentedMatrix();
                // Print Matrix
                printf("\nMatrix:\n");
                printMatrix(matrix);
                // GAUSS-SEIDAL
                applyGaussSeidal(matrix);
                // Free memory
                freeMatrix(matrix);
                break;
            case 7:
                root = getFunction(input, sizeof(input));
                // Print function tree
                printf("Function: ");
                printTree(root);
                printf("\n");
                // DIFFERENTIATION
                applyDifferentiation(root);
                // Free memory
                freeTree(root);
                break;
            case 8:
                root = getFunction(input, sizeof(input));
                // Print function tree
                printf("Function: ");
                printTree(root);
                printf("\n");
                // SIMPSON 1/3 AND 3/8
                applySimpson(root);
                // Free memory
                freeTree(root);
                break;
            case 9:
                root = getFunction(input, sizeof(input));
                // Print function tree
                printf("Function: ");
                printTree(root);
                printf("\n");
                // TRAPEZOIDAL
                rootValue = applyTrapezoidal(root);
                // Free memory
                freeTree(root);
                break;
            case 10:
                // GREGORY-NEWTON FORWARD INTERPOLATION
                applyGregoryNewtonForward();
                break;
        }
    }
    printf("\n\n");

    return 0;
}

/*
@brief Parses a mathematical expression string into an expression tree

@param str Pointer to the string containing the mathematical expression

@return Pointer to the root node of the parsed expression tree
*/
Node* parseFunction(const char **str) {
    int j;
    // Skip whitespace
    while (**str == ' ' || **str == '\t') (*str)++;

    // Read function name (max 9 characters)
    char func_name[10] = {0};
    const char *start = *str;
    int i = 0;
    while (isalpha(**str) && i < 9) {
        func_name[i++] = **str;
        (*str)++;
    }
    func_name[i] = '\0';

    // Convert func_name to lowercase for case-insensitive comparison
    for (j = 0; func_name[j]; j++) {
        func_name[j] = tolower((unsigned char)func_name[j]);
    }

    Node *log_base_node = NULL;
    // Check for log function base
    if (strcmp(func_name, "log") == 0 && **str == '_') {
        (*str)++; // Skip '_' character
        // Parse base expression
        if (**str == '(') {
            (*str)++; // Skip '('
            log_base_node = parseExpression(str);
            if (**str == ')') (*str)++; // Skip ')'
        } else {
            // If no parentheses, check simple cases first
            const char *base_start = *str;
            int is_simple = 1;
            int has_multiple_chars = 0;
            
            // Check if it's a simple case
            while (**str && **str != '(') {
                if ((**str < '0' || **str > '9') && 
                    (**str < 'a' || **str > 'z') && 
                    (**str < 'A' || **str > 'Z')) {
                    is_simple = 0;
                }
                if (**str != *base_start) {
                    has_multiple_chars = 1;
                }
                (*str)++;
            }
            
            *str = base_start;  // Reset pointer to start
            
            if (is_simple && !has_multiple_chars) {
                // Simple case - single character (number or variable)
                if (**str >= '0' && **str <= '9') {
                    Node *base_node = malloc(sizeof(Node));
                    base_node->type = NODE_NUMBER;
                    base_node->number = strtod(*str, (char**)str);
                    log_base_node = base_node;
                } else if ((**str >= 'a' && **str <= 'z') || (**str >= 'A' && **str <= 'Z')) {
                    if (**str == 'e' || **str == 'E') {
                        // Special case for e
                        Node *base_node = malloc(sizeof(Node));
                        base_node->type = NODE_NUMBER;
                        base_node->number = E;
                        (*str)++;
                        log_base_node = base_node;
                    } else {
                        Node *base_node = malloc(sizeof(Node));
                        base_node->type = NODE_VARIABLE;
                        base_node->variable = **str;
                        (*str)++;
                        log_base_node = base_node;
                    }
                }
            } else {
                // Complex case - parse entire expression
                const char *base_end = *str;
                char *base_str = malloc(base_end - base_start + 1);
                if (base_str) {
                    strncpy(base_str, base_start, base_end - base_start);
                    base_str[base_end - base_start] = '\0';
                    const char *temp_str = base_str;
                    log_base_node = parseExpression(&temp_str);
                    free(base_str);
                }
            }
        }
    }

    // If function name is found and is a supported function
    if (i > 0 && (
        strcmp(func_name, "sin") == 0 ||
        strcmp(func_name, "cos") == 0 ||
        strcmp(func_name, "tan") == 0 ||
        strcmp(func_name, "cot") == 0 ||
        strcmp(func_name, "cosec") == 0 ||
        strcmp(func_name, "sec") == 0 ||
        strcmp(func_name, "arcsin") == 0 ||
        strcmp(func_name, "arccos") == 0 ||
        strcmp(func_name, "arctan") == 0 ||
        strcmp(func_name, "arccot") == 0 ||
        strcmp(func_name, "arccosec") == 0 ||
        strcmp(func_name, "arcsec") == 0 ||
        strcmp(func_name, "log") == 0
    )) {
        Node *node = malloc(sizeof(Node));
        if (!node) return NULL;  // Check for memory error
        
        node->type = NODE_FUNCTION;
        strncpy(node->func.func, func_name, 9);
        node->func.func[9] = '\0';
        node->func.base = 0;  // Numeric base no longer used
        node->func.var_base = '\0';  // Variable base no longer used
        node->func.base_node = log_base_node;  // New field for base expression
        
        // Skip parentheses if present
        while (**str == ' ' || **str == '\t') (*str)++;
        if (**str == '(') {
            (*str)++; // Skip '('
            node->func.child = parseExpression(str);
            if (**str == ')') (*str)++; // Skip ')'
        } else {
            // If no parentheses, read single factor
            node->func.child = parseFactor(str);
        }
        return node;
    } else {
        // If no function name or not supported, reset pointer and parse as factor
        *str = start;
        return parseFactor(str);
    }
}

/*
@brief Parses a factor in a mathematical expression (numbers, variables, parentheses)

@param str Pointer to the string containing the factor

@return Pointer to the root node of the parsed factor
*/
Node* parseFactor(const char **str) {
    // Skip whitespace
    while (**str == ' ' || **str == '\t') (*str)++;
    
    Node *result = NULL;
    
    if (**str >= '0' && **str <= '9') {
        // Parse number
        Node *node = malloc(sizeof(Node));
        node->type = NODE_NUMBER;
        node->number = strtod(*str, (char**)str);
        result = node;
    } else if (**str == '(') {
        (*str)++; // Skip '('
        Node *node = parseExpression(str);
        if (**str == ')') (*str)++;
        result = node;
    } else if ((**str >= 'a' && **str <= 'z') || (**str >= 'A' && **str <= 'Z')) {
        // Parse variable or constant
        if (**str == 'e') {
            Node *node = malloc(sizeof(Node));
            node->type = NODE_NUMBER;
            node->number = E; // Value of e
            (*str)++;
            result = node;
        } else if (**str == 'p' && (*(*str + 1) == 'i' || *(*str + 1) == 'I')) {
            Node *node = malloc(sizeof(Node));
            node->type = NODE_NUMBER;
            node->number = PI; // Value of pi
            (*str) += 2; // Skip both 'p' and 'i'
            result = node;
        } else {
            Node *node = malloc(sizeof(Node));
            node->type = NODE_VARIABLE;
            node->variable = **str;
            (*str)++;
            result = node;
        }
    }

    // Check for implicit multiplication: if a number or parenthesis is followed by a variable
    if (result && (**str >= 'a' && **str <= 'z') || (**str >= 'A' && **str <= 'Z')) {
        Node *right = parseFactor(str);  // Parse the right variable
        if (right) {
            Node *mult = malloc(sizeof(Node));
            mult->type = NODE_OPERATOR;
            mult->op.op = '*';
            mult->op.left = result;
            mult->op.right = right;
            result = mult;
        }
    }
    
    return result;
}

/*
@brief Parses a power expression (exponentiation)

@param str Pointer to the string containing the power expression

@return Pointer to the root node of the parsed power expression
*/
Node* parsePower(const char **str) {
    Node *left = parseFunction(str);
    
    // Skip whitespace
    while (**str == ' ' || **str == '\t') (*str)++;
    
    while (**str == '^') {
        char op = **str;
        (*str)++;
        
        // Skip whitespace
        while (**str == ' ' || **str == '\t') (*str)++;
        
        Node *node = malloc(sizeof(Node));
        node->type = NODE_OPERATOR;
        node->op.op = op;
        node->op.left = left;
        // Exponentiation should bind from right, so we call recursively
        node->op.right = parsePower(str);
        left = node;
    }
    return left;
}

/*
@brief Parses a term in a mathematical expression (multiplication and division)

@param str Pointer to the string containing the term

@return Pointer to the root node of the parsed term
*/
Node* parseTerm(const char **str) {
    Node *left = parsePower(str);
    
    while (left && (**str == '*' || **str == '/' || 
           ((**str >= 'a' && **str <= 'z') || (**str >= 'A' && **str <= 'Z')) ||  // Variable
           (**str >= '0' && **str <= '9') ||  // Number
           **str == '(' ||  // Parenthesis
           isalpha(**str))) {  // Function
        
        // Skip whitespace
        while (**str == ' ' || **str == '\t') (*str)++;
        
        if (**str == '*' || **str == '/') {
            char op = **str;
            (*str)++;
            Node *node = malloc(sizeof(Node));
            node->type = NODE_OPERATOR;
            node->op.op = op;
            node->op.left = left;
            node->op.right = parsePower(str);
            left = node;
        } else {
            // Implicit multiplication or function
            Node *right = parsePower(str);
            if (right) {
                Node *node = malloc(sizeof(Node));
                node->type = NODE_OPERATOR;
                node->op.op = '*';
                node->op.left = left;
                node->op.right = right;
                left = node;
            }
        }
    }
    return left;
}

/*
@brief Parses a complete mathematical expression

@param str Pointer to the string containing the expression

@return Pointer to the root node of the parsed expression
*/
Node* parseExpression(const char **str) {
    Node *left = parseTerm(str);
    
    while (**str == '+' || **str == '-') {
        char op = **str;
        (*str)++;
        
        // Skip whitespace
        while (**str == ' ' || **str == '\t') (*str)++;
        
        Node *right = parseTerm(str);
        if (right) {
            Node *node = malloc(sizeof(Node));
            node->type = NODE_OPERATOR;
            node->op.op = op;
            node->op.left = left;
            node->op.right = right;
            left = node;
        }
    }
    return left;
}

/*
@brief Prints the expression tree in a readable format

@param node Pointer to the root node of the expression tree
*/
void printTree(Node *node) { // inorder
    if (!node) return;
    switch (node->type) {
        case NODE_NUMBER:
            printf("%g", node->number);
            break;
        case NODE_VARIABLE:
            printf("%c", node->variable);
            break;
        case NODE_OPERATOR:
            printf("(");
            printTree(node->op.left);
            printf(" %c ", node->op.op);
            printTree(node->op.right);
            printf(")");
            break;
        case NODE_FUNCTION:
            if (strcmp(node->func.func, "log") == 0) {
                printf("log_");
                if (node->func.base > 0) {
                    printf("%d", node->func.base);
                } else if (node->func.var_base != '\0') {
                    printf("%c", node->func.var_base);
                } else if (node->func.base_node != NULL) {
                    printTree(node->func.base_node);
                }
                printf("(");
                printTree(node->func.child);
                printf(")");
            } else {
                printf("%s(", node->func.func);
                printTree(node->func.child);
                printf(")");
            }
            break;
    }
}

/*
@brief Frees all memory allocated for the expression tree

@param node Pointer to the root node of the expression tree
*/
void freeTree(Node *node) {
    if (!node) return;
    switch (node->type) {
        case NODE_OPERATOR:
            freeTree(node->op.left);
            freeTree(node->op.right);
            break;
        case NODE_FUNCTION:
            freeTree(node->func.child);
            freeTree(node->func.base_node);
            break;
        default:
            break;
    }
    free(node);
}

/*
@brief Evaluates the expression tree for a given value of x

@param node Pointer to the root node of the expression tree
@param x The value to substitute for variables in the expression

@return The numerical result of evaluating the expression
*/
double evaluate(Node *node, double x) {
    if (!node) return 0.0;
    
    switch (node->type) {
        case NODE_NUMBER:
            return node->number;
        case NODE_VARIABLE:
            return x;  // Değişken yerine x değerini kullan
        case NODE_OPERATOR:
            switch (node->op.op) {
                case '+':
                    return evaluate(node->op.left, x) + evaluate(node->op.right, x);
                case '-':
                    return evaluate(node->op.left, x) - evaluate(node->op.right, x);
                case '*':
                    return evaluate(node->op.left, x) * evaluate(node->op.right, x);
                case '/': {
                    double denom = evaluate(node->op.right, x);
                    if (denom == 0.0) {
                        printf("Error: Division by zero.\n");
                        return 0.0;
                    }
                    return evaluate(node->op.left, x) / denom;
                }
                case '^': {
                    double base = evaluate(node->op.left, x);
                    double exp = evaluate(node->op.right, x);

                    // Negative base and fractional exponent → invalid
                    if (base < 0 && fmod(exp, 1.0) != 0.0) {
                        printf("Error: Negative number to fractional power is not defined in real numbers.\n");
                        return 0.0;
                    }

                    return pow(base, exp);
                }
                default:
                    return 0.0;
            }
        case NODE_FUNCTION: {
            double arg = evaluate(node->func.child, x);
            if (strcmp(node->func.func, "sin") == 0) {
                return sin(arg);
            } else if (strcmp(node->func.func, "cos") == 0) {
                return cos(arg);
            } else if (strcmp(node->func.func, "tan") == 0) {
                if (cos(arg) == 0.0) {
                    printf("Error: tan(%g) is not defined.\n", arg);
                    return 0.0;
                }
                return tan(arg);
            } else if (strcmp(node->func.func, "cot") == 0) {
                if (sin(arg) == 0.0) {
                    printf("Error: cot(%g) is not defined.\n", arg);
                    return 0.0;
                }
                return 1.0 / tan(arg);
            } else if (strcmp(node->func.func, "cosec") == 0) {
                if (sin(arg) == 0.0) {
                    printf("Error: cosec(%g) is not defined.\n", arg);
                    return 0.0;
                }
                return 1.0 / sin(arg);
            } else if (strcmp(node->func.func, "sec") == 0) {
                if (cos(arg) == 0.0) {
                    printf("Error: sec(%g) is not defined.\n", arg);
                    return 0.0;
                }
                return 1.0 / cos(arg);
            } else if (strcmp(node->func.func, "arcsin") == 0) {
                if (arg < -1.0 || arg > 1.0) {
                    printf("Error: arcsin(%g) is not defined (range: [-1,1]).\n", arg);
                    return 0.0;
                }
                return asin(arg);
            } else if (strcmp(node->func.func, "arccos") == 0) {
                if (arg < -1.0 || arg > 1.0) {
                    printf("Error: arccos(%g) is not defined (range: [-1,1]).\n", arg);
                    return 0.0;
                }
                return acos(arg);
            } else if (strcmp(node->func.func, "arctan") == 0) {
                return atan(arg);
            } else if (strcmp(node->func.func, "arccot") == 0) {
                return atan(1.0 / arg);
            } else if (strcmp(node->func.func, "arccosec") == 0) {
                if (arg == 0.0 || fabs(arg) < 1.0) {
                    printf("Error: arccosec(%g) is not defined (|x| >= 1).\n", arg);
                    return 0.0;
                }
                return asin(1.0 / arg);
            } else if (strcmp(node->func.func, "arcsec") == 0) {
                if (arg == 0.0 || fabs(arg) < 1.0) {
                    printf("Error: arcsec(%g) is not defined (|x| >= 1).\n", arg);
                    return 0.0;
                }
                return acos(1.0 / arg);
            } else if (strcmp(node->func.func, "log") == 0) {
                double base;
                if (node->func.base_node != NULL) {
                    base = evaluate(node->func.base_node, x);
                } else {
                    base = E;  // Default to natural logarithm
                }
                if (base <= 0.0 || base == 1.0) {
                    printf("Error: log base must be positive and different from 1.\n");
                    return 0.0;
                }
                if (arg <= 0.0) {
                    printf("Error: log(%g) is not defined (x > 0).\n", arg);
                    return 0.0;
                }
                return log(arg) / log(base);
            }
            return 0.0;
        }
    }
    return 0.0;
}

/*
@brief Displays the main menu and gets user's choice

@param choice Pointer to store the user's menu choice
*/
void printMenu(int* choice) {
    printf("0. Exit\n");
    printf("1. Bisection\n");
    printf("2. Regula-Falsi\n");
    printf("3. Newton-Raphson\n");
    printf("4. Inverse of Square Matrix\n");
    printf("5. Cholesky (ALU)\n");
    printf("6. Gauss Seidal\n");
    printf("7. Numerical Differentiation\n");
    printf("8. Simpson (1/3 and 3/8)\n");
    printf("9. Trapez\n");
    printf("10. Gregory-Newton Enterpolation\n");
    
    printf("Your choice: ");
    scanf("%d", choice);
    printf("\n");
}

/*
@brief Gets a mathematical function from user input

@param input Buffer to store the input string
@param size Size of the input buffer

@return Pointer to the root node of the parsed function
*/
Node* getFunction(char* input, size_t size) {
    printf("Enter function: ");
    // Clear input buffer
    int c;
    while ((c = getchar()) != '\n' && c != EOF);
    // Get function input
    fgets(input, size, stdin);
    const char *p = input;
    Node *root = parseExpression(&p);
    return root;
}

/*
@brief Compares signs of two numbers

@param a First number
@param b Second number

@return 1 if signs are same, 0 if different
*/
int compareSigns(double a, double b) {
    if ((a >= 0 && b >= 0) || (a <= 0 && b <= 0)) {
        return 1;
    }
    return 0;
}

/*
@brief Finds a root of a function using the bisection method

@param root Root node of the function tree
@param start Start of the interval
@param end End of the interval
@param epsilon Desired accuracy
@param maxIterations Maximum number of iterations

@return The approximate root value
*/
double bisection(Node* root, double start, double end, double epsilon, int maxIterations) {
    double mid;
    int iterations = 0;
    int epsilonCond = 1;
    double error;
    double midResult, startResult, endResult;

    startResult = evaluate(root, start);
    endResult = evaluate(root, end);

    printf("\nBISECTION\n");

    if (compareSigns(startResult, endResult) == 1) {
        printf("There is no real root in the given range of the function.\n\n");
        return 0;
    }

    do {
        iterations++;
        error = (end - start); // / pow(2, iterations);
        if (error < epsilon) {
            epsilonCond = 0;
        }

        mid = (start + end) / 2;
        midResult = evaluate(root, mid);

        printf("Iteration %d -> Range: [%.4f, %.4f], x%d = %.3f, f(x%d) = %.3f\n", iterations, start, end, iterations, mid, iterations, midResult);

        // Check if we found a root (function value is exactly zero)
        if (midResult == 0.0) {
            printf("Root found: x = %.4f (function value is exactly zero)\n\n", mid);
            return mid;
        }

        if (compareSigns(startResult, midResult) == 1) {
            start = mid;
            startResult = midResult;
        }
        else {
            end = mid;
            endResult = midResult;
        }
    }
    while (iterations < maxIterations && epsilonCond);

    if (epsilonCond == 0) {
        printf("Root found: x = %.4f (error is less than epsilon)\n\n", mid);
    }
    else {
        printf("Maximum number of iterations was reached. Root found: x = %.4f\n\n", mid);
    }

    return mid;
}

/*
@brief Applies the bisection method with user input

@param root Root node of the function tree
@return The approximate root value
*/
double applyBisection(Node* root) {
    double start, end, epsilon;
    int maxIterations;

    printf("Enter start value: ");
    scanf("%lf", &start);
    printf("Enter end value: ");
    scanf("%lf", &end);
    printf("Enter epsilon value: ");
    scanf("%lf", &epsilon);
    printf("Enter max iterations: ");
    scanf("%d", &maxIterations);

    return bisection(root, start, end, epsilon, maxIterations);
}

/*
@brief Finds a root of a function using the regula falsi method

@param root Root node of the function tree
@param start Start of the interval
@param end End of the interval
@param epsilon Desired accuracy
@param maxIterations Maximum number of iterations

@return The approximate root value
*/
double regulaFalsi(Node* root, double start, double end, double epsilon, int maxIterations) {
    double c;
    int iterations = 0;
    int epsilonCond = 1;
    double error;
    double cResult, startResult, endResult;

    startResult = evaluate(root, start);
    endResult = evaluate(root, end);

    printf("\nREGULA-FALSI\n");

    if (compareSigns(startResult, endResult) == 1) {
        printf("There is no real root in the given range of the function.\n\n");
        return 0;
    }

    do {
        iterations++;
        error = (end - start); // / pow(2, iterations);
        if (error < epsilon) {
            epsilonCond = 0;
        }

        c = (start * endResult - end * startResult) / (endResult - startResult);
        cResult = evaluate(root, c);

        // Check if we found a root (function value is exactly zero)
        if (cResult == 0.0) {
            printf("Root found: x = %.4f (function value is exactly zero)\n\n", c);
            return c;
        }

        printf("Iteration %d -> Range: [%.4f, %.4f], x%d = %.3f, f(x%d) = %.3f\n", iterations, start, end, iterations, c, iterations, cResult);

        // Check if we found a root (function value is exactly zero)
        if (cResult == 0.0) {
            printf("Root found: x = %.4f (function value is exactly zero)\n\n", c);
            return c;
        }

        if (compareSigns(startResult, cResult) == 1) {
            start = c;
            startResult = cResult;
        }
        else {
            end = c;
            endResult = cResult;
        }
    }
    while (iterations < maxIterations && epsilonCond);

    if (epsilonCond == 0) {
        printf("Root found: x = %.4f (error is less than epsilon)\n\n", c);
    }
    else {
        printf("Maximum number of iterations was reached. Root found: x = %.4f\n\n", c);
    }

    return c;
}

/*
@brief Applies the regula falsi method with user input

@param root Root node of the function tree

@return The approximate root value
*/
double applyRegulaFalsi(Node* root) {
    double start, end, epsilon;
    int maxIterations;

    printf("Enter start value: ");
    scanf("%lf", &start);
    printf("Enter end value: ");
    scanf("%lf", &end);
    printf("Enter epsilon value: ");
    scanf("%lf", &epsilon);
    printf("Enter max iterations: ");
    scanf("%d", &maxIterations);

    return regulaFalsi(root, start, end, epsilon, maxIterations);
}

/*
@brief Calculates the derivative using forward difference method

@param root Root node of the function tree
@param x Point at which to calculate derivative
@param h Step size

@return The approximate derivative value
*/
double forwardDifferenceDerivative(Node* root, double x, double h) {
    double f1 = evaluate(root, x + h);
    double f2 = evaluate(root, x);
    return (f1 - f2) / h;
}

/*
@brief Calculates the derivative using backward difference method

@param root Root node of the function tree
@param x Point at which to calculate derivative
@param h Step size

@return The approximate derivative value
*/
double backwardDifferenceDerivative(Node* root, double x, double h) {
    double f1 = evaluate(root, x);
    double f2 = evaluate(root, x - h);
    return (f1 - f2) / h;
}

/*
@brief Calculates the derivative using central difference method

@param root Root node of the function tree
@param x Point at which to calculate derivative
@param h Step size

@return The approximate derivative value
*/
double centralDifferenceDerivative(Node* root, double x, double h) {
    double f1 = evaluate(root, x + h);
    double f2 = evaluate(root, x - h);
    return (f1 - f2) / (2 * h);
}

/*
@brief Applies numerical differentiation with user input

@param root Root node of the function tree
*/
void applyDifferentiation(Node* root) {
    double x, h;

    printf("Enter the point of differentiation: ");
    scanf("%lf", &x);
    printf("Enter the step size: ");
    scanf("%lf", &h);

    printf("\nForward Difference Derivative:\n");
    printf("f'(%.4f) = (%.3f - %.3f) / %.3f = %.3f\n", x, evaluate(root, x + h), evaluate(root, x), h, forwardDifferenceDerivative(root, x, h));

    printf("Backward Difference Derivative:\n");
    printf("f'(%.4f) = (%.3f - %.3f) / %.3f = %.3f\n", x, evaluate(root, x), evaluate(root, x - h), h, backwardDifferenceDerivative(root, x, h));

    printf("Central Difference Derivative:\n");
    printf("f'(%.4f) = (%.3f - %.3f) / (2 * %.3f) = %.3f\n", x, evaluate(root, x + h), evaluate(root, x - h), h, centralDifferenceDerivative(root, x, h));

    printf("\n");
}

/*
@brief Finds a root of a function using the Newton-Raphson method

@param root Root node of the function tree
@param start Start of the search range
@param end End of the search range
@param x0 Initial guess
@param epsilon Desired accuracy
@param maxIterations Maximum number of iterations

@return The approximate root value
*/
double newtonRaphson(Node* root, double x0, double epsilon, int maxIterations) {
    double currentX, nextX;
    int iterations = 0;
    int epsilonCond = 1;
    double start = x0 - 10, end = x0 + 10;
    double error, currentXResult, currentDerivative;
    double prevError = INFINITY;  // Previous error
    int divergenceCount = 0;      // Divergence counter
    const int MAX_DIVERGENCE = 3; // Maximum divergence attempts
    const double DIVERGENCE_THRESHOLD = 1.5; // Divergence threshold

    currentX = x0;
    currentXResult = evaluate(root, currentX);
    currentDerivative = centralDifferenceDerivative(root, currentX, 0.0001);

    printf("\nNEWTON-RAPHSON\n");

    do {
        iterations++;

        nextX = currentX - (currentXResult / currentDerivative);

        error = fabs(nextX - currentX);
        
        // Divergence check
        if (error > prevError * DIVERGENCE_THRESHOLD) {
            divergenceCount++;
            if (divergenceCount >= MAX_DIVERGENCE) {
                // Random new starting point
                double newX0 = start + ((double)rand() / RAND_MAX) * (end - start);
                printf("Divergence detected. New starting point: %.4f\n", newX0);
                return newtonRaphson(root, newX0, epsilon, maxIterations);
            }
        }

        if (error < epsilon) {
            epsilonCond = 0;
        }

        printf("Iteration %d -> x%d = %.3f - (%.3f / %.3f) = %.3f\n", 
               iterations, iterations, currentX, currentXResult, currentDerivative, nextX);

        currentX = nextX;
        currentXResult = evaluate(root, currentX);
        currentDerivative = centralDifferenceDerivative(root, currentX, 0.0001);
        prevError = error;
    }
    while (iterations < maxIterations && epsilonCond);

    if (epsilonCond == 0) {
        printf("Root found: x = %.4f (error is less than epsilon)\n\n", nextX);
    }
    else {
        printf("Maximum number of iterations was reached. Root found: x = %.4f\n\n", nextX);
    }

    return nextX;
}

/*
@brief Applies the Newton-Raphson method with user input

@param root Root node of the function tree

@return The approximate root value
*/
double applyNewtonRaphson(Node* root) {
    double x0, epsilon;
    int maxIterations;

    printf("Enter the starting point: ");
    scanf("%lf", &x0);   
    printf("Enter epsilon value: ");
    scanf("%lf", &epsilon);
    printf("Enter max iterations: ");
    scanf("%d", &maxIterations);

    return newtonRaphson(root, x0, epsilon, maxIterations);
}

/*
@brief Performs numerical integration using Simpson's 1/3 rule

@param root Root node of the function tree
@param start Start of the interval
@param end End of the interval
@param n Number of intervals

@return The approximate integral value
*/
double simpson1_3(Node* root, double start, double end, double n) {
    double h, currentX, startValue, endValue;
    double sum = 0;
    int iteration = 0;

    h = (end - start) / n;
    currentX = start + h;

    printf("\nSIMPSON 1/3\n");

    startValue = evaluate(root, start);
    endValue = evaluate(root, end);

    sum += startValue + endValue;

    printf("x%d - > f(%.4f) = %.4f\n", iteration, start, startValue);
    while (currentX != end) {
        iteration++;
        printf("x%d - > f(%.4f) = %.4f\n", iteration, currentX, evaluate(root, currentX));
        if (iteration % 2 == 0) {
            sum += 2 * evaluate(root, currentX);
        }
        else {
            sum += 4 * evaluate(root, currentX);
        }
        currentX += h;
    }
    printf("x%d - > f(%.4f) = %.4f\n", iteration + 1, end, endValue);

    printf("Result: (%.2lf / 3) * (%.4lf) = %.4lf\n\n", h, sum, (h / 3) * sum);
    return (h / 3) * sum;
}

/*
@brief Applies Simpson's 1/3 rule with user input

@param root Root node of the function tree
@return The approximate integral value
*/
double applySimpson1_3(Node* root) {
    double start, end, n;

    printf("Enter the starting point: ");
    scanf("%lf", &start);
    printf("Enter the ending point: ");
    scanf("%lf", &end);
    printf("Enter the number of intervals: ");
    scanf("%lf", &n);
    
    return simpson1_3(root, start, end, n);
}

/*
@brief Performs numerical integration using Simpson's 3/8 rule

@param root Root node of the function tree
@param start Start of the interval
@param end End of the interval
@param n Number of intervals

@return The approximate integral value
*/
double simpson3_8(Node* root, double start, double end, double n) {
    double limitsH = (end - start) / n;
    double h = limitsH / 3;
    int iteration = 1;
    double newStart, newEnd;
    double result, sum = 0;
    int j, loopIteration = 0;
    double values[4];
    double* results;

    results = (double*)calloc(n, sizeof(double));

    printf("\nSIMPSON 3/8\n");

    for (iteration = 1; iteration <= n; iteration++) {
        sum = 0;
        printf("%d. Iteration\n", iteration);

        newStart = start + ((iteration - 1) * limitsH);
        newEnd = newStart + limitsH;

        for (j = 0; j <= 3; j++) {
            printf("x%d -> f(%.4f) = %.4f\n", j, newStart + (j * limitsH / 3), evaluate(root, newStart + (j * limitsH / 3)));
            values[j] = evaluate(root, newStart + (j * limitsH / 3));
            if (j == 0 || j == 3) {
                sum += values[j];
            }
            else {
                sum += 3 * values[j];
            }
        }

        results[iteration - 1] = (limitsH / 8) * sum;  
        printf("Result: (%.4f - %.4f) * (%.4f + 3 * %.4f + 3 * %.4f + %.4f) / 8 = %.4f\n\n", newEnd, newStart, values[0], values[1], values[2], values[3], results[iteration - 1]);
        printf("\n");
    }

    printf("Total Result = ");
    for (j = 0; j < n; j++) {
        if (j < n - 1) {
            printf("%.4f + ", results[j]);
        }
        else {
            printf("%.4f = ", results[j]);
        }
        result += results[j];
    }
    printf("%.4f\n\n", result);

    return result;
}

/*
@brief Applies Simpson's 3/8 rule with user input

@param root Root node of the function tree

@return The approximate integral value
*/
double applySimpson3_8(Node* root) {
    double start, end, n;

    printf("Enter the starting point: ");
    scanf("%lf", &start);
    printf("Enter the ending point: ");
    scanf("%lf", &end);
    printf("Enter the number of intervals (n): ");
    scanf("%lf", &n);

    return simpson3_8(root, start, end, n);
}

/*
@brief Applies both Simpson's 1/3 and 3/8 rules with user input

@param root Root node of the function tree
@return Array containing results from both methods
*/
double* applySimpson(Node* root) {
    double start, end;
    double simpson38_n;
    double simpson13_n;
    int n = 2;
    double* results = calloc(n, sizeof(double));
    
    if (results == NULL) {
    	printf("Error in memory allocation!");
    	return NULL;
	}

    printf("Enter the starting point: ");
    scanf("%lf", &start);
    printf("Enter the ending point: ");
    scanf("%lf", &end);

    // Simpson 1/3
    printf("Enter the number of intervals for Simpson 1/3: ");
    scanf("%lf", &simpson13_n);
    results[0] = simpson1_3(root, start, end, simpson13_n);

    // Simpson 3/8
    printf("Enter the number of intervals for Simpson 3/8: ");
    scanf("%lf", &simpson38_n);
    results[1] = simpson3_8(root, start, end, simpson38_n);

    return results;
}

/*
@brief Performs numerical integration using the trapezoidal rule

@param root Root node of the function tree
@param start Start of the interval
@param end End of the interval
@param n Number of intervals

@return The approximate integral value
*/
double trapezoidal(Node* root, double start, double end, double n) {
    double h = (end - start) / n;
    double sum = 0;
    int i;
    double* values;

    values = (double*)calloc(n + 1, sizeof(double));

    printf("\nTRAPEZOIDAL\n");

    for (i = 0; i <= n; i++) {
        values[i] = evaluate(root, start + i * h);
        printf("x%d -> f(%.4f) = %.4f\n", i, start + i * h, values[i]);
        if (i != n) {
            sum += 2 * values[i];
        }
        else {
            sum += values[i];
        }
    }

    printf("Result: (%.4f / 2) * (", h);
    for (i = 1; i < n + 1; i++) {
        if (i == 1) {
            printf("%.4f + ", values[i]);
        }
        else if (i == n) {
            printf("%.4f)", values[i]);
        }
        else {
            printf("2 * %.4f + ", values[i]);
        }
    }
    printf(") = %.4f\n\n", (h / 2) * sum);

    free(values);

    return (h / 2) * sum;
}

/*
@brief Applies the trapezoidal rule with user input

@param root Root node of the function tree

@return The approximate integral value
*/
double applyTrapezoidal(Node* root) {
    double start, end, n;

    printf("Enter the starting point: ");
    scanf("%lf", &start);
    printf("Enter the ending point: ");
    scanf("%lf", &end);
    printf("Enter the number of intervals (n): ");
    scanf("%lf", &n);

    return trapezoidal(root, start, end, n);
}

/*
@brief Creates a new matrix with specified dimensions

@param rows Number of rows
@param cols Number of columns

@return The created matrix
*/
Matrix createMatrix(int rows, int cols) {
    Matrix m;
    int i;

    m.rows = rows;
    m.cols = cols;

    m.data = (double**)calloc(rows, sizeof(double *));
    for(i = 0; i < rows; i++) {
        m.data[i] = (double*)calloc(cols, sizeof(double));
    }

    return m;
}

/*
@brief Frees memory allocated for a matrix

@param m The matrix to free
*/
void freeMatrix(Matrix m) {
    int i;

    for(i = 0; i < m.rows; i++) {
        free(m.data[i]);
    }
    free(m.data);
}

/*
@brief Gets matrix dimensions and elements from user input

@return The created matrix
*/
Matrix getMatrix() {
    Matrix m;
    int rows, cols, i, j;

    printf("Enter the number of rows: ");
    scanf("%d", &rows);
    printf("Enter the number of columns: ");
    scanf("%d", &cols);

    m = createMatrix(rows, cols);
    for(i = 0; i < rows; i++) {
        for(j = 0; j < cols; j++) {
            printf("Enter the element at position [%d][%d]: ", i + 1, j + 1);
            scanf("%lf", &m.data[i][j]);
        }
    }

    return m;
}

/*
@brief Gets an augmented matrix from user input

@return The created augmented matrix
*/
Matrix getAugmentedMatrix() {
    Matrix m;
    int rows, cols, i, j;

    printf("Enter the number of rows of your augmented matrix: ");
    scanf("%d", &rows);

    cols = rows + 1;
    m = createMatrix(rows, cols);
    for(i = 0; i < rows; i++) {
        for(j = 0; j < cols; j++) {
            printf("Enter the element at position [%d][%d]: ", i + 1, j + 1);
            scanf("%lf", &m.data[i][j]);
        }
    }

    return m;
}

/*
@brief Gets a square matrix from user input

@return The created square matrix
*/
Matrix getSquareMatrix() {
    Matrix m;
    int rows, cols, i, j;

    printf("Enter the number of rows of your square matrix: ");
    scanf("%d", &rows);

    cols = rows;
    m = createMatrix(rows, cols);
    for(i = 0; i < rows; i++) {
        for(j = 0; j < cols; j++) {
            printf("Enter the element at position [%d][%d]: ", i + 1, j + 1);
            scanf("%lf", &m.data[i][j]);
        }
    }

    return m;
}

/*
@brief Multiplies two matrices

@param a First matrix
@param b Second matrix

@return The product matrix
*/
Matrix multiplyMatrices(Matrix a, Matrix b) {
    if (a.cols != b.rows) {
        printf("Matrix dimensions do not match for multiplication.\n");
        return (Matrix){NULL, 0, 0};
    }

    Matrix result = createMatrix(a.rows, b.cols);
    int i, j, k;

    for(i = 0; i < a.rows; i++) {
        for(j = 0; j < b.cols; j++) {
            result.data[i][j] = 0;
            for(k = 0; k < a.cols; k++) {
                result.data[i][j] += a.data[i][k] * b.data[k][j];
            }
        }
    }

    return result;
}

/*
@brief Separates an augmented matrix into coefficient and result matrices

@param a The augmented matrix

@return Array containing coefficient and result matrices
*/
Matrix* separateAugmentedMatrix(Matrix a) {
    Matrix coefs, results;
    int i, j;
    Matrix* result;

    coefs = createMatrix(a.rows, a.cols - 1);
    results = createMatrix(a.rows, 1);

    for (i = 0; i < a.rows; i++) {
        for (j = 0; j < a.cols - 1; j++) {
            coefs.data[i][j] = a.data[i][j];
        }
        results.data[i][0] = a.data[i][a.cols - 1];
    }

    result = (Matrix*)malloc(sizeof(Matrix) * 2);
    result[0] = coefs;
    result[1] = results;

    return result;
}

/*
@brief Checks if a matrix is square

@param a The matrix to check

@return 1 if square, 0 otherwise
*/
int isSquareMatrix(Matrix a) {
    return a.rows == a.cols;
}

/*
@brief Checks if a matrix is augmented

@param a The matrix to check

@return 1 if augmented, 0 otherwise
*/
int isAugmentedMatrix(Matrix a) {
    return a.rows + 1 == a.cols;
}

/*
@brief Prints an augmented matrix equation

@param a The augmented matrix
@param title Title for the equation
@param var Variable name

@param indexed Whether to use indexed variables
*/
void printAugmentedMatrixEquation(Matrix a, char* title, char var, int indexed) {
    Matrix* matrices;
    int i, j;
    char varName[3];  // enough for x0, x1, x2 or a, b, c

    if (!isAugmentedMatrix(a)) {
        printf("Matrix is not augmented.\n");
        return;
    }

    matrices = separateAugmentedMatrix(a);
    
    printf("\n%s:\n", title);
    printf("----------------\n");
    
    // Coefficient matrix
    for(i = 0; i < matrices[0].rows; i++) {
        printf("| ");
        for(j = 0; j < matrices[0].cols; j++) {
            printf("%8.4f ", matrices[0].data[i][j]);
        }
        printf("|");
        
        // Variables matrix
        printf("   | ");
        if(i < matrices[0].cols) {
            if(indexed) {
                printf("%c%d ", tolower(var), i);
            } else {
                printf("%c ", var + (i % 26));
            }
        } else {
            printf("  ");
        }
        printf("|");
        
        // Equal sign
        if(i == matrices[0].rows / 2) {
            printf(" = ");
        } else {
            printf("   ");
        }
        
        // Results matrix
        printf("| %8.4f |\n", matrices[1].data[i][0]);
    }
    printf("----------------\n\n");

    freeMatrix(matrices[0]);
    freeMatrix(matrices[1]);
    free(matrices);
}

/*
@brief Prints a matrix equation

@param coefs Coefficient matrix
@param results Result matrix
@param title Title for the equation
@param var Variable name

@param indexed Whether to use indexed variables
*/
void printMatrixEquation(Matrix coefs, Matrix results, char* title, char var, int indexed) {
    int i, j;
    char varName[3];  // enough for x0, x1, x2 or a, b, c

    if (!isSquareMatrix(coefs)) {
        printf("Matrix is not square.\n");
        return;
    }
    
    printf("\n%s:\n", title);
    printf("----------------\n");
    
    // Coefficient matrix
    for(i = 0; i < coefs.rows; i++) {
        printf("| ");
        for(j = 0; j < coefs.cols; j++) {
            printf("%8.4f ", coefs.data[i][j]);
        }
        printf("|");
        
        // Variables matrix
        printf("   | ");
        if(i < coefs.cols) {
            if(indexed) {
                printf("%c%d ", tolower(var), i);
            } else {
                printf("%c ", var + (i % 26));
            }
        } else {
            printf("  ");
        }
        printf("|");
        
        // Equal sign
        if(i == coefs.rows / 2) {
            printf(" = ");
        } else {
            printf("   ");
        }
        
        // Results matrix
        printf("| %8.4f |\n", results.data[i][0]);
    }
    printf("----------------\n\n");
}

/*
@brief Performs Gaussian elimination to convert a matrix to upper triangular form

@param a The input matrix
@param det Pointer to store the determinant value
@param rank Pointer to store the rank of the matrix

@return The upper triangular form of the input matrix
*/
Matrix upperTriangularGaussianElimination(Matrix a, double* det, int* rank) {
    int n, i, k, j, maxRow, signChanges = 0;
    double temp, factor, pivot;
    Matrix upper;

    n = a.rows;
    *det = 1.0;
    *rank = n;
    
    // Create a copy of the matrix
    upper = createMatrix(n, n);
    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            upper.data[i][j] = a.data[i][j];
        }
    }

    // Gaussian elimination
    for(i = 0; i < n; i++) {
        // Pivot selection
        maxRow = i;
        for(k = i + 1; k < n; k++) {
            if(fabs(upper.data[k][i]) > fabs(upper.data[maxRow][i])) {
                maxRow = k;
            }
        }

        // Row swapping
        if(maxRow != i) {
            for(j = 0; j < n; j++) {
                double temp = upper.data[i][j];
                upper.data[i][j] = upper.data[maxRow][j];
                upper.data[maxRow][j] = temp;
            }
            signChanges++;
        }

        // If pivot element is 0, rank decreases
        if(fabs(upper.data[i][i]) < 1e-10) {
            (*rank)--;
            continue;
        }

        // Normalize pivot row
        pivot = upper.data[i][i];
        *det *= pivot;
        for(j = i; j < n; j++) {
            upper.data[i][j] /= pivot;
        }

        // Elimination
        for(k = i + 1; k < n; k++) {
            factor = upper.data[k][i];
            for(j = i; j < n; j++) {
                upper.data[k][j] -= factor * upper.data[i][j];
            }
        }
    }

    // Determinant sign
    *det *= (signChanges % 2 == 0) ? 1 : -1;

    return upper;
}

/*
@brief Performs Gaussian elimination to convert a matrix to lower triangular form

@param a The input matrix
@param det Pointer to store the determinant value
@param rank Pointer to store the rank of the matrix

@return The lower triangular form of the input matrix
*/
Matrix lowerTriangularGaussianElimination(Matrix a, double* det, int* rank) {
    int n, i, k, j, maxRow, signChanges = 0;
    double temp, factor, pivot;
    Matrix lower;

    n = a.rows;
    *det = 1.0;
    *rank = n;
    
    // Create a copy of the matrix
    lower = createMatrix(n, n);
    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            lower.data[i][j] = a.data[i][j];
        }
    }

    // Gaussian elimination (for lower triangular form)
    for(i = n-1; i >= 0; i--) {
        // Pivot selection
        maxRow = i;
        for(k = i - 1; k >= 0; k--) {
            if(fabs(lower.data[k][i]) > fabs(lower.data[maxRow][i])) {
                maxRow = k;
            }
        }

        // Row swapping
        if(maxRow != i) {
            for(j = 0; j < n; j++) {
                temp = lower.data[i][j];
                lower.data[i][j] = lower.data[maxRow][j];
                lower.data[maxRow][j] = temp;
            }
            signChanges++;
        }

        // If pivot element is zero, decrease rank
        if(fabs(lower.data[i][i]) < 1e-10) {
            (*rank)--;
            continue;
        }

        // Normalize pivot row (make diagonal 1)
        pivot = lower.data[i][i];
        *det *= pivot;
        for(j = 0; j <= i; j++) {
            lower.data[i][j] /= pivot;
        }

        // Elimination (subtract from rows above)
        for(k = i - 1; k >= 0; k--) {
            factor = lower.data[k][i];
            for(j = 0; j <= i; j++) {
                lower.data[k][j] -= factor * lower.data[i][j];
            }
        }
    }

    // Determinant sign
    *det *= (signChanges % 2 == 0) ? 1 : -1;

    return lower;
}

/*
@brief Solves a system of linear equations using an augmented matrix

@param a The augmented matrix

@return Array containing the solution values
*/
double* solveAugmentedMatrix(Matrix a) {
    Matrix* separated;
    Matrix coefs, results;
    int i, j, rank;
    double det, sum;
    double* variables;

    // Separate augmented matrix into coefficient and result matrices
    separated = separateAugmentedMatrix(a);
    coefs = separated[0];
    results = separated[1];

    // Check if coefficient matrix is square
    if (!isSquareMatrix(coefs)) {
        printf("Coefficient matrix must be square.\n");
        free(separated);
        return NULL;
    }

    // Allocate memory for variables
    variables = (double*)calloc(coefs.rows, sizeof(double));
    if (!variables) {
        printf("Memory allocation failed.\n");
        free(separated);
        return NULL;
    }

    // Convert to lower triangular form
    coefs = lowerTriangularGaussianElimination(coefs, &det, &rank);
    
    // Check rank
    if (rank < coefs.rows) {
        printf("System is underdetermined or inconsistent.\n");
        free(variables);
        free(separated);
        return NULL;
    }

    // Calculate first variable
    variables[0] = results.data[0][0];

    // Calculate other variables (for lower triangular form)
    for (i = 1; i < coefs.rows; i++) {
        sum = 0.0;
        for (j = 0; j < i; j++) {
            sum += coefs.data[i][j] * variables[j];
        }
        variables[i] = results.data[i][0] - sum;
    }

    // Clean up memory
    freeMatrix(coefs);
    freeMatrix(results);
    free(separated);

    return variables;
}

/*
@brief Solves a system of linear equations

@param coefs Coefficient matrix
@param results Result matrix

@return Array containing the solution values
*/
double* solveEquationSystem(Matrix coefs, Matrix results) {
    int i, j, rank = coefs.rows;  // Initialize rank to number of rows
    double det, sum;
    double* variables;

    // Check if coefficient matrix is square
    if (!isSquareMatrix(coefs)) {
        printf("Coefficient matrix must be square.\n");
        return NULL;
    }

    // Allocate memory for variables
    variables = (double*)calloc(coefs.rows, sizeof(double));
    if (!variables) {
        printf("Memory allocation failed.\n");
        return NULL;
    }

    // Check if matrix is upper or lower triangular
    int isUpper = 1;
    int isLower = 1;
    
    for(i = 0; i < coefs.rows; i++) {
        for(j = 0; j < coefs.cols; j++) {
            if(i > j && fabs(coefs.data[i][j]) > 1e-10) isUpper = 0;
            if(i < j && fabs(coefs.data[i][j]) > 1e-10) isLower = 0;
        }
    }

    // If neither upper nor lower triangular, convert to lower triangular
    if(!isUpper && !isLower) {
        coefs = lowerTriangularGaussianElimination(coefs, &det, &rank);
        isLower = 1;
    }

    // Check rank
    for(i = 0; i < coefs.rows; i++) {
        if(fabs(coefs.data[i][i]) < 1e-10) {
            rank--;
        }
    }
    
    if (rank < coefs.rows) {
        printf("System is underdetermined or inconsistent.\n");
        free(variables);
        return NULL;
    }

    if(isUpper) {
        // Back substitution for upper triangular
        for(i = coefs.rows - 1; i >= 0; i--) {
            sum = 0.0;
            for(j = i + 1; j < coefs.cols; j++) {
                sum += coefs.data[i][j] * variables[j];
            }
            variables[i] = (results.data[i][0] - sum) / coefs.data[i][i];
        }
    } else {
        // Forward substitution for lower triangular
        variables[0] = results.data[0][0] / coefs.data[0][0];
        for(i = 1; i < coefs.rows; i++) {
            sum = 0.0;
            for(j = 0; j < i; j++) {
                sum += coefs.data[i][j] * variables[j];
            }
            variables[i] = (results.data[i][0] - sum) / coefs.data[i][i];
        }
    }

    return variables;
}

/*
@brief Checks if a matrix is symmetric

@param a The matrix to check

@return 1 if symmetric, 0 otherwise
*/
int isSymmetric(Matrix a) {
    if (!isSquareMatrix(a)) return 0;
    
    int i, j;
    for(i = 0; i < a.rows; i++) {
        for(j = 0; j < i; j++) {
            if(fabs(a.data[i][j] - a.data[j][i]) > 1e-10) {
                return 0;
            }
        }
    }
    return 1;
}

/*
@brief Checks if a matrix is positive definite

@param a The matrix to check

@return 1 if positive definite, 0 otherwise
*/
int isPositiveDefinite(Matrix a) {
    if (!isSymmetric(a)) return 0;
    
    int n = a.rows;
    int i, j, k;
    double sum;
    Matrix temp = createMatrix(n, n);
    
    // Matrisin bir kopyasını oluştur
    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            temp.data[i][j] = a.data[i][j];
        }
    }
    
    // Cholesky ayrıştırması sırasında kontrol
    for(i = 0; i < n; i++) {
        // Köşegen elemanı için kontrol
        sum = temp.data[i][i];
        for(k = 0; k < i; k++) {
            sum -= temp.data[i][k] * temp.data[i][k];
        }
        if(sum <= 1e-10) {
            freeMatrix(temp);
            return 0;
        }
        temp.data[i][i] = sqrt(sum);
        
        // Diğer elemanlar için kontrol
        for(j = i + 1; j < n; j++) {
            sum = temp.data[j][i];
            for(k = 0; k < i; k++) {
                sum -= temp.data[i][k] * temp.data[j][k];
            }
            temp.data[j][i] = sum / temp.data[i][i];
        }
    }
    
    freeMatrix(temp);
    return 1;
}

/*
@brief Calculates the sum of numbers from 1 to n

@param n Upper limit

@return The sum
*/
double sumUpTo(int n) {
    return n * (n + 1) / 2;
}

/*
@brief Calculates the factorial of a number

@param n The number
@return The factorial value
*/
double factorial(int n) {
    int i;
    double result = 1;
    for (i = 1; i <= n; i++) {
        result *= i;
    }
    return result;
}

/*
@brief Calculates the transpose of a matrix

@param a The input matrix

@return The transposed matrix
*/
Matrix getTranspose(Matrix a) {
    int i, j;
    Matrix transpose;

    transpose = createMatrix(a.cols, a.rows);

    for (i = 0; i < a.rows; i++) {
        for(j = 0; j < a.cols; j++) {
            transpose.data[j][i] = a.data[i][j];
        }
    }

    return transpose;
}

/*
@brief Performs Cholesky factorization on a positive definite matrix

@param a The input matrix

@return Array containing the lower and upper triangular matrices
*/
Matrix* choleskyFactorization(Matrix a) {
    int size;
    Matrix* matrices;
    int i, j, k;
    double sum;

    if (!isPositiveDefinite(a)) {
        printf("Matrix is not positive definite.\n\n");
        return NULL;
    }

    size = sumUpTo(a.rows);
    matrices = (Matrix*)malloc(2 * sizeof(Matrix));

    // Lower ve Upper matrislerini oluştur
    matrices[0] = createMatrix(a.rows, a.rows); // Lower
    matrices[1] = createMatrix(a.rows, a.rows); // Upper

    // Cholesky ayrıştırması
    for(i = 0; i < a.rows; i++) {
        for(j = 0; j <= i; j++) {
            sum = a.data[i][j];
            
            // Köşegen elemanlar için
            if(i == j) {
                for(k = 0; k < j; k++) {
                    sum -= matrices[0].data[i][k] * matrices[0].data[i][k];
                }
                matrices[0].data[i][j] = sqrt(sum);
            }
            // Köşegen altı elemanlar için
            else {
                for(k = 0; k < j; k++) {
                    sum -= matrices[0].data[i][k] * matrices[0].data[j][k];
                }
                matrices[0].data[i][j] = sum / matrices[0].data[j][j];
            }
        }
    }

    // Upper (Transpose of Lower)
    matrices[1] = getTranspose(matrices[0]);

    return matrices;
}

/*
@brief Prints matrix decomposition in a formatted way

@param title Title for the decomposition
@param A Original matrix
@param a Name for matrix A
@param L Lower triangular matrix
@param l Name for matrix L
@param U Upper triangular matrix
@param u Name for matrix U
*/
void printMatrixDecomposition(char* title, Matrix A, char* a, Matrix L, char* l, Matrix U, char* u) {
    int i, j;
    
    printf("\n%s:\n", title);
    
    printf("----------------\n");
    // Her satır için
    for(i = 0; i < A.rows; i++) {
        // A matrisi
        printf("| ");
        for(j = 0; j < A.cols; j++) {
            printf("%8.4f ", A.data[i][j]);
        }
        printf("|");
        
        // İlk satırda = işareti
        if(i == 0) {
            printf(" = ");
        } else {
            printf("   ");
        }
        
        // L matrisi
        printf("| ");
        for(j = 0; j < L.cols; j++) {
            printf("%8.4f ", L.data[i][j]);
        }
        printf("|");
        
        // İlk satırda * işareti
        if(i == 0) {
            printf(" * ");
        } else {
            printf("   ");
        }
        
        // U matrisi
        printf("| ");
        for(j = 0; j < U.cols; j++) {
            printf("%8.4f ", U.data[i][j]);
        }
        printf("|\n");
    }
    printf("----------------\n");
}

/*
@brief Performs Cholesky decomposition on a positive definite matrix

@param a The input matrix

@return Array containing the solution values
*/
double* choleskyDecomposition(Matrix a) {
    Matrix* matrices;
    Matrix* LU;
    Matrix newResults;
    double* results;
    int i, j, k;
    double sum;

    // Separate augmented matrix into coefficient and result matrices
    matrices = separateAugmentedMatrix(a);
    
    // Check if the matrix is positive definite
    if (!isPositiveDefinite(matrices[0])) {
        printf("Error: Matrix is not positive definite. Cholesky decomposition cannot be performed.\n\n");
        freeMatrix(matrices[0]);
        freeMatrix(matrices[1]);
        free(matrices);
        return NULL;
    }
    
    // Perform Cholesky factorization
    LU = choleskyFactorization(matrices[0]);
    results = (double*)malloc(a.rows * sizeof(double));

    // Print matrix decomposition
    printMatrixDecomposition("A = L * U Decomposition", matrices[0], "A", LU[0], "L", LU[1], "U");
    
    // Solve L * y = b
    printMatrixEquation(LU[0], matrices[1], "L = Results Matrix", 'a', 0);
    results = solveEquationSystem(LU[0], matrices[1]);
    for (i = 0; i < a.rows; i++) {
        printf("%c = %lf\n", 97 + i, results[i]);
    }

    // Create new results matrix for U * x = y
    newResults = createMatrix(a.rows, 1);
    for (i = 0; i < a.rows; i++) {
        newResults.data[i][0] = results[i];
    }

    // Solve U * x = y
    printMatrixEquation(LU[1], newResults, "U = L's Results Matrix", 'x', 1);
    results = solveEquationSystem(LU[1], newResults);

    // Clean up memory
    freeMatrix(LU[0]);
    freeMatrix(LU[1]);
    free(LU);
    freeMatrix(matrices[0]);
    freeMatrix(matrices[1]);
    free(matrices);
    freeMatrix(newResults);

    return results;
}

/*
@brief Applies Cholesky decomposition with user input

@param a The input matrix
*/
void applyCholesky(Matrix a) {
    double* results;
    int i;

    results = choleskyDecomposition(a);
    if (results != NULL) {
        printf("Results:\n");
        for (i = 0; i < a.rows; i++) {
            printf("x%d = %lf\n", i, results[i]);
        }
        printf("\n");
        free(results);
    }
}

/*
@brief Calculates the determinant of a square matrix

@param a The input matrix

@return The determinant value
*/
double determinant(Matrix a) {
    double det;
    int rank;
    
    if (!isSquareMatrix(a)) {
        printf("Matrix is not square.\n");
        return 0;
    }

    Matrix upper = upperTriangularGaussianElimination(a, &det, &rank);

    freeMatrix(upper);

    return det;
}

/*
@brief Checks if a matrix is invertible

@param a The matrix to check

@return 1 if invertible, 0 otherwise
*/
int isInvertible(Matrix a) {
    return a.rows == a.cols && determinant(a) != 0;
}

/*
@brief Calculates the inverse of a square matrix

@param a The input matrix

@return The inverse matrix
*/
Matrix inverseMatrix(Matrix a) {
    Matrix inverse, aug;
    int rank, maxRow, i, j, k, n = a.rows;
    double det, temp, pivot, factor;

    if (!isSquareMatrix(a)) {
        printf("Matrix is not square.\n");
        return (Matrix){NULL, 0, 0};
    }
    
    // Create augmented matrix [A|I]
    aug = createMatrix(n, 2*n);
    
    // Copy original matrix to left side
    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            aug.data[i][j] = a.data[i][j];
        }
    }
    
    // Add identity matrix to right side
    for(i = 0; i < n; i++) {
        for(j = n; j < 2*n; j++) {
            aug.data[i][j] = (i == j-n) ? 1.0 : 0.0;
        }
    }

    // Gauss-Jordan elimination
    for(i = 0; i < n; i++) {
        // Pivot selection
        maxRow = i;
        for(k = i + 1; k < n; k++) {
            if(fabs(aug.data[k][i]) > fabs(aug.data[maxRow][i])) {
                maxRow = k;
            }
        }

        // If pivot element is 0, matrix is not invertible
        if(fabs(aug.data[maxRow][i]) < 1e-10) {
            printf("Matrix is not invertible.\n");
            freeMatrix(aug);
            return (Matrix){NULL, 0, 0};
        }

        // Row swapping
        if(maxRow != i) {
            for(j = 0; j < 2*n; j++) {
                temp = aug.data[i][j];
                aug.data[i][j] = aug.data[maxRow][j];
                aug.data[maxRow][j] = temp;
            }
        }

        // Normalize pivot row
        pivot = aug.data[i][i];
        for(j = 0; j < 2*n; j++) {
            aug.data[i][j] /= pivot;
        }

        // Subtract pivot row from other rows
        for(k = 0; k < n; k++) {
            if(k != i) {
                factor = aug.data[k][i];
                for(j = 0; j < 2*n; j++) {
                    aug.data[k][j] -= factor * aug.data[i][j];
                }
            }
        }
    }

    // Separate inverse matrix
    inverse = createMatrix(n, n);
    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            inverse.data[i][j] = aug.data[i][j+n];
        }
    }

    freeMatrix(aug);

    return inverse;
}

/*
@brief Prints a matrix in a formatted way

@param m The matrix to print
*/
void printMatrix(Matrix m) {
    int i, j;
    for(i = 0; i < m.rows; i++) {
        for(j = 0; j < m.cols; j++) {
            printf("%8.4f ", m.data[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

/*
@brief Multiplies a matrix with its transpose

@param a The input matrix

@return The product matrix
*/
Matrix multiplyMatrixWithTranspose(Matrix a) {
    int n = a.rows;
    Matrix result = createMatrix(n, n);
    int i, j, k;
    
    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            result.data[i][j] = 0;
            for(k = 0; k < n; k++) {
                result.data[i][j] += a.data[i][k] * a.data[j][k];
            }
        }
    }
    
    return result;
}

/*
@brief Swaps two rows in a matrix

@param a Pointer to the matrix
@param row1 First row index
@param row2 Second row index
*/
void swapRows(Matrix* a, int row1, int row2) {
    double* temp;
    temp = a->data[row1];
    a->data[row1] = a->data[row2];
    a->data[row2] = temp;
}

/*
@brief Solves a system of linear equations using the Gauss-Seidel method

@param a The augmented matrix
@param epsilon Desired accuracy
@param maxIterations Maximum number of iterations

@return Array containing the solution values
*/
double* gaussSeidal(Matrix a, double epsilon, int maxIterations) {
    Matrix* matrices;
    Matrix coefs, results;
    int i, j, k, iterations = 0;
    int epsilonCond = 0;
    double max, sum, oldValue;
    double *values, *oldValues;

    // Separate augmented matrix
    matrices = separateAugmentedMatrix(a);
    coefs = matrices[0];
    results = matrices[1];

    // Allocate memory for current and previous values
    values = (double*)calloc(coefs.cols, sizeof(double));
    oldValues = (double*)calloc(coefs.cols, sizeof(double));

    // Rearrange rows for dominant diagonal
    for (i = 0; i < coefs.rows; i++) {
        max = -DBL_MAX;
        for (j = 0; j < coefs.cols; j++) {
            if (fabs(coefs.data[i][j]) > max) {
                max = fabs(coefs.data[i][j]);
                k = j;
            }
        }
        swapRows(&coefs, i, k);
        swapRows(&results, i, k);
    }

    printf("\n");
    printMatrix(coefs);

    // Get initial values from user
    for(i = 0; i < coefs.cols; i++) {
        printf("Enter initial value for x%d: ", i);
        scanf("%lf", &values[i]);
    }
    printf("\n");

    // Print first iteration results
    printf("Iteration %d:\n", iterations + 1);
    for(i = 0; i < coefs.cols; i++) {
        printf("x%d = %lf, deltax%d = %lf\n", 
                i, values[i], i, fabs(values[i] - oldValues[i]));
    }
    printf("\n");
    iterations++;

    // Main iteration loop
    while (!epsilonCond && iterations < maxIterations) {
        // Save previous values
        for(i = 0; i < coefs.cols; i++) {
            oldValues[i] = values[i];
        }

        // Update each variable
        for(i = 0; i < coefs.rows; i++) {
            sum = 0.0;
            for(j = 0; j < coefs.cols; j++) {
                if(j != i) {
                    sum += coefs.data[i][j] * values[j];
                }
            }
            values[i] = (results.data[i][0] - sum) / coefs.data[i][i];
        }

        // Check convergence
        epsilonCond = 1;
        for(i = 0; i < coefs.cols; i++) {
            if(fabs(values[i] - oldValues[i]) > epsilon) {
                epsilonCond = 0;
                break;
            }
        }

        // Print iteration results
        printf("Iteration %d:\n", iterations + 1);
        for(i = 0; i < coefs.cols; i++) {
            printf("x%d = %lf, deltax%d = %lf\n", 
                   i, values[i], i, fabs(values[i] - oldValues[i]));
        }
        printf("\n");

        iterations++;
    }

    // Print final results
    printf("Results:\n");
    for(i = 0; i < coefs.cols; i++) {
        printf("x%d = %lf\n", i, values[i]);
    }
    printf("\n");

    // Clean up memory
    free(oldValues);
    freeMatrix(matrices[0]);
    freeMatrix(matrices[1]);
    free(matrices);

    return values;
}

/*
@brief Applies the Gauss-Seidel method with user input

@param matrix The input matrix

@return Array containing the solution values
*/
double* applyGaussSeidal(Matrix matrix) {
    double epsilon;
    int maxIterations;
    double* results;

    printf("Enter epsilon: ");
    scanf("%lf", &epsilon);
    printf("Enter max iterations: ");
    scanf("%d", &maxIterations);

    return gaussSeidal(matrix, epsilon, maxIterations);
}

/*
@brief Creates a forward difference table for interpolation

@param x Array of x values
@param y Array of y values
@param n Number of points

@return 2D array containing the difference table
*/
double** getForwardDifferenceTable(double* x, double* y, int n) {
    int i, j;
    double h = x[1] - x[0];  // Step size
    double** deltas;  // 2D array for difference table
    
    // Allocate memory for difference table
    deltas = (double**)malloc(n * sizeof(double*));
    for(i = 0; i < n; i++) {
        deltas[i] = (double*)calloc(n, sizeof(double));
    }
    
    // Fill first column with y values
    for(i = 0; i < n; i++) {
        deltas[i][0] = y[i];
    }
    
    // Calculate forward differences
    for(j = 1; j < n; j++) {
        for(i = 0; i < n-j; i++) {
            deltas[i][j] = deltas[i+1][j-1] - deltas[i][j-1];
        }
    }
    
    return deltas;
}

/*
@brief Prints a forward difference table

@param x Array of x values
@param arr The difference table
@param n Number of points
*/
void printForwardDifferenceTable(double* x, double** arr, int n) {
    int i, j;
    char header[20];
    int nonZeroCols = 0;

    // Count non-zero columns
    for(j = 1; j < n; j++) {
        if(fabs(arr[0][j]) > 1e-10) {
            nonZeroCols++;
        }
    }

    // Print table header
    printf("\nForward Difference Table:\n");
    printf("%-12s %-12s ", "x", "y");
    for(i = 1; i <= nonZeroCols; i++) {
        sprintf(header, "D^%dy", i);
        printf("%-12s ", header);
    }
    printf("\n");
    
    // Print separator line
    for(i = 0; i <= nonZeroCols + 1; i++) {
        printf("------------");
    }
    printf("\n");
    
    // Print table contents
    for(i = 0; i < n; i++) {
        printf("%-12.4f ", x[i]);
        for(j = 0; j < n-i; j++) {
            if(j == 0 || fabs(arr[i][j]) > 1e-10) {
                printf("%-12.4f ", arr[i][j]);
            }
        }
        printf("\n");
    }
    
    // Print bottom separator
    for(i = 0; i <= nonZeroCols + 1; i++) {
        printf("------------");
    }
    printf("\n");
}

/*
@brief Calculates the product of differences for interpolation

@param x Array of x values
@param size Size of the array
@param searchX The x value to interpolate
@param iteration Current iteration

@return The product value
*/
double product(double* x, int size, double searchX, int iteration) {
    int i;
    double product = 1;

    // Calculate product of differences
    for (i = 0; i < iteration; i++) {
        if (i < size) {
            product *= (searchX - x[i]);
        }       
    }

    return product;
}

/*
@brief Prints the product terms for interpolation

@param x Array of x values
@param size Size of the array
@param searchX The x value to interpolate
@param iteration Current iteration
*/
void printProducts(double* x, int size, double searchX, int iteration) {
    int i;

    // Print product terms
    for (i = 0; i < iteration; i++) {
        if (i == iteration - 1) {
            printf("(%.3lf - %.3lf)", searchX, x[i]);
        }
        else {
            printf("(%.3lf - %.3lf) * ", searchX, x[i]);
        }    
    }
}

/*
@brief Performs Gregory-Newton forward interpolation

@param x Array of x values
@param y Array of y values
@param n Number of points
@param searchX The x value to interpolate

@return The interpolated value
*/
double gregoryNewtonForward(double* x, double* y, int n, double searchX) {
    int i, j;
    double h, sum;
    double** deltas;

    // Calculate step size
    h = x[1] - x[0];
    
    // Generate and print difference table
    deltas = getForwardDifferenceTable(x, y, n);
    printForwardDifferenceTable(x, deltas, n);

    // Print interpolation formula
    printf("\nGREGORY-NEWTON INTERPOLATION\n");
    printf("F(X) = %.3lf", y[0]);
    sum = y[0];
    
    // Calculate interpolated value
    for (i = 1; i < n; i++) {
        if(fabs(deltas[0][i]) > 1e-10) {
            printf(" + (%.3lf * ", deltas[0][i]);
            printProducts(x, n, searchX, i);
            printf(") / (%.3lf^%d * %d!)", h, i, i);
            sum += (deltas[0][i] * product(x, n, searchX, i)) / (pow(h, i) * factorial(i));
        }
    }
    printf("\n");

    // Print final result
    printf("F(%.3lf) = %.3lf\n\n", searchX, sum);

    // Clean up memory
    for (i = 0; i < n; i++) {
        free(deltas[i]);
    }
    free(deltas);
    
    return sum;
}

/*
@brief Applies Gregory-Newton forward interpolation with user input

@return The interpolated value
*/
double applyGregoryNewtonForward() {
    double *x, *y;
    int i, n;
    double x0, searchX, h;

    // Get number of points
    printf("Enter the number of points: ");
    scanf("%d", &n);

    // Allocate memory for x and y arrays
    x = (double*)malloc(n * sizeof(double));
    y = (double*)malloc(n * sizeof(double));

    // Get initial x point and step size
    printf("Enter the first x point: ");
    scanf("%lf", &x0);
    printf("Enter the step size (h): ");
    scanf("%lf", &h);

    // Generate x values
    for (i = 0; i < n; i++) {
        x[i] = x0 + (i * h);
    }

    // Get y values
    for (i = 0; i < n; i++) {
        printf("Enter y%d: ", i);
        scanf("%lf", &y[i]);
    }

    // Get interpolation point
    printf("Enter the x point you want to learn it's value: ");
    scanf("%lf", &searchX);

    // Perform interpolation
    double result = gregoryNewtonForward(x, y, n, searchX);

    // Clean up memory
    free(x);
    free(y);

    return result;
}
