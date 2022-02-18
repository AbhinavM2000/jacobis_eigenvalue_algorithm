#include<stdio.h>

#include<math.h>


void printline(int k) //Function to print a line
{
  int i;
  for (i = 1; i <= k; i++) {
    printf("_");
  }
  printf("\n");
}

void readmatrix(float * x, int a) //Function to read a square matrix of order 'a'
{
  int i, j;
  printf("Enter %d elements = ", a * a);
  for (i = 1; i <= a; i++) {
    for (j = 1; j <= a; j++) {
      scanf("%f", & x[i * 10 + j]);
    }
  }
}

void matrixequals(float * x, float * y, int n) { //Function to do x=y for matrices

  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= n; j++) {
      x[i * 10 + j] = y[i * 10 + j];
    }
  }

}

void writematrix(float * x, int a) //Function to print an square matrix of order 'a'
{
  int i, j;
  for (i = 1; i <= a; i++) {
    for (j = 1; j <= a; j++) {
      printf("%f\t", x[i * 10 + j]);
    }
    printf("\n");
  }
}

void transpose(float * a, int n, float * trans) //Function to transpose an square matrix of order 'n'
{
  int i, j;
  for (i = 1; i <= n; i++) {
    for (j = 1; j <= n; j++) {
      trans[i * 10 + j] = a[j * 10 + i];
    }
  }
}

void mulply(float * x, float * y, float * res, int n) //Function for matrix multiplication
{
  int i, j, k;
  for (i = 1; i <= n; i++) {
    for (j = 1; j <= n; j++) {
      res[i * 10 + j] = 0;
      for (k = 1; k <= n; k++) {
        res[i * 10 + j] = res[i * 10 + j] + x[i * 10 + k] * y[k * 10 + j];
      }
    }
  }

}

void rotmat(float C, float S, int n, int l, int k, float * a) //Function to construct a rotation matrix from Jacobi Method's algorithm
{
  int i, j;
  for (i = 1; i <= n; i++) {
    for (j = 1; j <= n; j++) {
      if (i == j) {
        a[i * 10 + j] = 1;
      }
      if (i != j) {
        a[i * 10 + j] = 0;
      }
    }
  }
  a[l * 10 + l] = C;
  a[l * 10 + k] = -S;
  a[k * 10 + l] = S;
  a[k * 10 + k] = C;
}

float J_AJ(float * x, float * y, int n, float * r) { //Function to compute J'A J 
  float x1[100], r1[100];
  transpose(x, n, x1);
  mulply(x1, y, r1, n);
  mulply(r1, x, r, n);
}

float clear_arr(float r[100], float b[100], int n, float a[100]) { //Function to clear values of all arrays
  int i, j;
  for (i = 1; i <= n; i++) {
    for (j = 1; j <= n; j++) {
      a[i * 10 + j] = r[i * 10 + j];
      r[i * 10 + j] = 0;
      b[i * 10 + j] = 0;
    }
  }
}

float sumod(float * z, int n) {   //Funtion to return the sum of squares of off-diagonal elements
  float sum = 0;
  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= n; j++) {
      if (i != j)
        sum = sum + fabs(z[i * 10 + j] * z[i * 10 + j]);
    }
  }

  return sum;
}

void main() {
  int i, j, k, n, p, q, num;
  num = 0;
  float a[100], b[100], r[100], r11[100], newb[100], multipliedpart[100], A, B, C, S, big, sum, M;
  printline(60);
  printf("Enter The order of the real symmetric matrix : ");
  scanf("%d", & n);
  printline(60);
  readmatrix(a, n);
  printline(60);
  printf("The given real symmetric matrix is :\n");
  writematrix(a, n);
  printline(60);
  printf("Enter the tolerance for sum of squares of off-diagonal elements (eg: 0.0001) : ");
  scanf("%f", & M);
  sum = sumod(a, n);
//Setting initial value of multipliedpart to identity matrix
  for (int ii = 1; ii <= n; ii++) {
    for (int jj = 1; jj <= n; jj++) {
      if (ii == jj) {
        multipliedpart[ii * 10 + jj] = 1;
      }
      if (ii != jj) {
        multipliedpart[ii * 10 + jj] = 0;
      }
    }
  }

  while (sum > M) {
    big = 0;
    for (i = 1; i <= n; i++) {
      for (j = 1; j <= n; j++) {
        if (i < j) {
          if (fabs(a[i * 10 + j]) > big) {
            big = fabs(a[i * 10 + j]);      //Finding the largest off-diagonal element
            p = i;
            q = j;
          }
        }
      }
    }

    num++;
    A = atan(a[p * 10 + q]) / (a[p * 10 + p] - a[q * 10 + q]); //angle of rotation
    C = cos(A);
    S = sin(A);
    rotmat(C, S, n, p, q, b);
//Conscecutively multiplying the J matrices to obtain eigenvalue matrix
    matrixequals(newb, b, n);
    mulply(multipliedpart, b, r11, n);
    matrixequals(multipliedpart, r11, n);
//Finding J'AJ
    J_AJ(b, a, n, r);
    sum = sumod(a, n);
    clear_arr(r, b, n, a);
  }
  printline(60);
  printf("The required eigen values are :\n");
  for (i = 1; i <= n; i++) {
    printf("a%d%d = %f.\n", i, i, a[i * 10 + i]);
  }

  printf("\nEigenvector matrix : \n");
  writematrix(multipliedpart, n);

  printf("\nNumer of rotations performed : %d\n", num);

  printline(60);
}
