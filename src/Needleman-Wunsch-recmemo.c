/**
 * \file Needleman-Wunsch-recmemo.c
 * \brief recursive implementation with memoization of Needleman-Wunsch global alignment algorithm that computes the distance between two genetic sequences
 * \version 0.1
 * \date 03/10/2022
 * \author Jean-Louis Roch (Ensimag, Grenoble-INP - University Grenoble-Alpes) jean-louis.roch@grenoble-inp.fr
 *
 * Documentation: see Needleman-Wunsch-recmemo.h
 * Costs of basic base opertaions (SUBSTITUTION_COST, SUBSTITUTION_UNKNOWN_COST, INSERTION_COST) are
 * defined in Needleman-Wunsch-recmemo.h
 */

#include "Needleman-Wunsch-recmemo.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* for strchr */
// #include <ctype.h> /* for toupper */
#include "characters_to_base.h" /* mapping from char to base */

/*****************************************************************************/

/* Context of the memoization : passed to all recursive calls */
/** \def NOT_YET_COMPUTED
 * \brief default value for memoization of minimal distance (defined as an impossible value for a distance, -1).
 */
#define NOT_YET_COMPUTED -1L

/** \struct NW_MemoContext
 * \brief data for memoization of recursive Needleman-Wunsch algorithm
 */
struct NW_MemoContext
{
   char *X;     /*!< the longest genetic sequences */
   char *Y;     /*!< the shortest genetic sequences */
   size_t M;    /*!< length of X */
   size_t N;    /*!< length of Y,  N <= M */
   long **memo; /*!< memoization table to store memo[0..M][0..N] (including stopping conditions phi(M,j) and phi(i,N) */
};

/*
 *  static long EditDistance_NW_RecMemo(struct NW_MemoContext *c, size_t i, size_t j)
 * \brief  EditDistance_NW_RecMemo :  Private (static)  recursive function with memoization \
 * direct implementation of Needleman-Wursch extended to manage FASTA sequences (cf TP description)
 * \param c : data passed for recursive calls that includes the memoization array
 * \param i : starting position of the left sequence :  c->X[ i .. c->M ]
 * \param j : starting position of the right sequence :  c->Y[ j .. c->N ]
 */
static long EditDistance_NW_RecMemo(struct NW_MemoContext *c, size_t i, size_t j)
/* compute and returns phi(i,j) using data in c -allocated and initialized by EditDistance_NW_Rec */
{
   if (c->memo[i][j] == NOT_YET_COMPUTED)
   {
      long res;
      char Xi = c->X[i];
      char Yj = c->Y[j];
      if (i == c->M) /* Reach end of X */
      {
         if (j == c->N)
            res = 0; /* Reach end of Y too */
         else
            res = (isBase(Yj) ? INSERTION_COST : 0) + EditDistance_NW_RecMemo(c, i, j + 1);
      }
      else if (j == c->N) /* Reach end of Y but not end of X */
      {
         res = (isBase(Xi) ? INSERTION_COST : 0) + EditDistance_NW_RecMemo(c, i + 1, j);
      }
      else if (!isBase(Xi)) /* skip ccharacter in Xi that is not a base */
      {
         ManageBaseError(Xi);
         res = EditDistance_NW_RecMemo(c, i + 1, j);
      }
      else if (!isBase(Yj)) /* skip ccharacter in Yj that is not a base */
      {
         ManageBaseError(Yj);
         res = EditDistance_NW_RecMemo(c, i, j + 1);
      }
      else
      {             /* Note that stopping conditions (i==M) and (j==N) are already stored in c->memo (cf EditDistance_NW_Rec) */
         long min = /* initialization  with cas 1*/
             (isUnknownBase(Xi) ? SUBSTITUTION_UNKNOWN_COST
                                : (isSameBase(Xi, Yj) ? 0 : SUBSTITUTION_COST)) +
             EditDistance_NW_RecMemo(c, i + 1, j + 1);
         {
            long cas2 = INSERTION_COST + EditDistance_NW_RecMemo(c, i + 1, j);
            if (cas2 < min)
               min = cas2;
         }
         {
            long cas3 = INSERTION_COST + EditDistance_NW_RecMemo(c, i, j + 1);
            if (cas3 < min)
               min = cas3;
         }
         res = min;
      }
      c->memo[i][j] = res;
   }
   return c->memo[i][j];
}

/* EditDistance_NW_Rec :  is the main function to call, cf .h for specification
 * It allocates and initailizes data (NW_MemoContext) for memoization and call the
 * recursivefunction EditDistance_NW_RecMemo
 * See .h file for documentation
 */
long EditDistance_NW_Rec(char *A, size_t lengthA, char *B, size_t lengthB)
{
   _init_base_match();
   struct NW_MemoContext ctx;
   if (lengthA >= lengthB) /* X is the longest sequence, Y the shortest */
   {
      ctx.X = A;
      ctx.M = lengthA;
      ctx.Y = B;
      ctx.N = lengthB;
   }
   else
   {
      ctx.X = B;
      ctx.M = lengthB;
      ctx.Y = A;
      ctx.N = lengthA;
   }
   size_t M = ctx.M;
   size_t N = ctx.N;
   { /* Allocation and initialization of ctx.memo to NOT_YET_COMPUTED*/
      /* Note: memo is of size (N+1)*(M+1) but is stored as (M+1) distinct arrays each with (N+1) continuous elements
       * It would have been possible to allocate only one big array memezone of (M+1)*(N+1) elements
       * and then memo as an array of (M+1) pointers, the memo[i] being the address of memzone[i*(N+1)].
       */
      ctx.memo = (long **)malloc((M + 1) * sizeof(long *));
      if (ctx.memo == NULL)
      {
         perror("EditDistance_NW_Rec: malloc of ctx_memo");
         exit(EXIT_FAILURE);
      }
      for (int i = 0; i <= M; ++i)
      {
         ctx.memo[i] = (long *)malloc((N + 1) * sizeof(long));
         if (ctx.memo[i] == NULL)
         {
            perror("EditDistance_NW_Rec: malloc of ctx_memo[i]");
            exit(EXIT_FAILURE);
         }
         for (int j = 0; j <= N; ++j)
            ctx.memo[i][j] = NOT_YET_COMPUTED;
      }
   }

   /* Compute phi(0,0) = ctx.memo[0][0] by calling the recursive function EditDistance_NW_RecMemo */
   long res = EditDistance_NW_RecMemo(&ctx, 0, 0);

   { /* Deallocation of ctx.memo */
      for (int i = 0; i <= M; ++i)
         free(ctx.memo[i]);
      free(ctx.memo);
   }
   return res;
}

void print_array(long *array, size_t length)
{
   printf("[");
   for (int i = 0; i < length; i++)
   {
      if (i == length - 1)
      {
         printf("%ld", array[i]);
      }
      else
      {

         printf("%ld, ", array[i]);
      }
   }
   printf("] len(%ld)\n", length);
}

long min(long a, long b, long c)
{
   if (a <= b && a <= c)
   {
      return a;
   }
   else if (b <= a && b <= c)
   {
      return b;
   }
   else
   {
      return c;
   }
}

long EditDistance_NW_Iter(char *A, size_t lengthA, char *B, size_t lengthB)
{
   _init_base_match();
   struct NW_MemoContext ctx;

   if (lengthA >= lengthB)
   {
      ctx.X = A;
      ctx.M = lengthA;
      ctx.Y = B;
      ctx.N = lengthB;
   }
   else
   {
      ctx.X = B;
      ctx.M = lengthB;
      ctx.Y = A;
      ctx.N = lengthA;
   }
   size_t M = ctx.M;
   size_t N = ctx.N;

   long *Y_col = (long *)malloc((N + 1) * sizeof(long));
   if (Y_col == NULL)
   {
      perror("EditDistance_NW_Iter: malloc of Y_col. You must have at least one!");
      exit(EXIT_FAILURE);
   }
   Y_col[N] = 0;
   for (int row = N - 1; row >= 0; row--)
   {
      Y_col[row] = (isBase(ctx.Y[row]) ? INSERTION_COST : 0) + Y_col[row + 1];
   }

   long prev_value;

   for (int col = M - 1; col >= 0; col--)
   {
      for (int row = N; row >= 0; row--)
      {
         // printf("row=%ld, col=%ld\n", row, col);
         if (row == N)
         {
            prev_value = Y_col[row];
            Y_col[row] = (isBase(ctx.X[col]) ? INSERTION_COST : 0) + Y_col[row];
         }
         else if (!isBase(ctx.X[col]))
         {
            prev_value = Y_col[row];
            ManageBaseError(ctx.X[col]);
         }
         else if (!isBase(ctx.Y[row]))
         {
            prev_value = Y_col[row];
            Y_col[row] = Y_col[row + 1];
            ManageBaseError(ctx.Y[row]);
         }
         else
         {
            long diag = (isUnknownBase(ctx.X[col]) ? SUBSTITUTION_UNKNOWN_COST : (isSameBase(ctx.X[col], ctx.Y[row]) ? 0 : SUBSTITUTION_COST)) + prev_value;
            long top = INSERTION_COST + Y_col[row + 1];
            long right = INSERTION_COST + Y_col[row];
            prev_value = Y_col[row];
            Y_col[row] = min(diag, right, top);
         }
         // print_array(Y_col, N + 1);
      }
   }
   free(Y_col);

   return Y_col[0];
}

long EditDistance_NW_Iter_CA(char *A, size_t lengthA, char *B, size_t lengthB)
{
   _init_base_match();
   struct NW_MemoContext ctx;

   if (lengthA >= lengthB)
   {
      ctx.X = A;
      ctx.M = lengthA;
      ctx.Y = B;
      ctx.N = lengthB;
   }
   else
   {
      ctx.X = B;
      ctx.M = lengthB;
      ctx.Y = A;
      ctx.N = lengthA;
   }
   size_t M = ctx.M;
   size_t N = ctx.N;

   long *Y_col = (long *)malloc((N + 1) * sizeof(long));
   if (Y_col == NULL)
   {
      perror("EditDistance_NW_Iter: malloc of Y_col. You must have at least one!");
      exit(EXIT_FAILURE);
   }
   Y_col[N] = 0;
   for (int row = N - 1; row >= 0; row--)
   {
      Y_col[row] = (isBase(ctx.Y[row]) ? INSERTION_COST : 0) + Y_col[row + 1];
   }
   long prev_value;
   long Z = 4096 / 2;
   long number_of_Z_rows = (N) / Z + 1;
   long number_of_Z_cols = (M) / Z + 1;
   long Z_row_length = Z > M ? M : Z;
   long *Z_row = (long *)malloc((Z_row_length + 1) * sizeof(long));
   long prev_z;

   for (int z_col = number_of_Z_cols - 1; z_col >= 0; z_col--)
   {
      long max_cols = M - (number_of_Z_cols - 1 - z_col) * Z > Z ? Z : M - (number_of_Z_cols - 1 - z_col) * Z;
      for (int z_row = number_of_Z_rows - 1; z_row >= 0; z_row--)
      {
         long max_rows = N - (number_of_Z_rows - 1 - z_row) * Z > Z ? Z : N - (number_of_Z_rows - 1 - z_row) * Z;
         if (z_row != number_of_Z_rows - 1)
         {
            max_rows--;
         }
         for (int counter_col = max_cols - 1; counter_col >= 0; counter_col--)
         {
            long col = z_col > 0 ? counter_col + (z_col - 1) * Z + M - (number_of_Z_cols - 1) * Z : counter_col;
            long value_to_keep;
            for (int counter_row = max_rows; counter_row >= 0; counter_row--)
            {
               long row = z_row > 0 ? counter_row + (z_row - 1) * Z + N - (number_of_Z_rows - 1) * Z : counter_row;

               if (row == N)
               {
                  prev_value = Y_col[row];
                  Y_col[row] = (isBase(ctx.X[col]) ? INSERTION_COST : 0) + Y_col[row];
               }
               else if (!isBase(ctx.X[col]))
               {
                  prev_value = Y_col[row];
                  ManageBaseError(ctx.X[col]);
               }
               else if (!isBase(ctx.Y[row]))
               {
                  prev_value = Y_col[row];
                  if (counter_row == max_rows)
                  {
                     Y_col[row] = Z_row[counter_col + Z_row_length - max_cols];
                  }
                  else
                  {

                     Y_col[row] = Y_col[row + 1];
                  }
                  ManageBaseError(ctx.Y[row]);
               }
               else
               {
                  if (z_row == number_of_Z_rows - 1 || counter_row != max_rows)
                  {

                     long diag = (isUnknownBase(ctx.X[col]) ? SUBSTITUTION_UNKNOWN_COST : (isSameBase(ctx.X[col], ctx.Y[row]) ? 0 : SUBSTITUTION_COST)) + prev_value;
                     long top = INSERTION_COST + Y_col[row + 1];
                     long right = INSERTION_COST + Y_col[row];
                     // printf("row=%ld, col=%ld, z_row=%ld, z_col=%ld, counter_row=%ld, counter_col=%ld\n", row, col, z_row, z_col, counter_row, counter_col);
                     // printf("prev_value=%ld\n", prev_value);
                     // printf("NORMAL diag=%ld, top=%ld, right=%ld\n", diag, top, right);
                     prev_value = Y_col[row];
                     Y_col[row] = min(diag, right, top);
                     // print_array(Y_col, N + 1);
                  }
                  else
                  {
                     // print_array(Z_row, Z_row_length + 1);
                     // print_array(Y_col, N + 1);
                     long diag;
                     if (counter_col == max_cols - 1)
                     {

                        diag = (isUnknownBase(ctx.X[col]) ? SUBSTITUTION_UNKNOWN_COST : (isSameBase(ctx.X[col], ctx.Y[row]) ? 0 : SUBSTITUTION_COST)) + Z_row[counter_col + Z_row_length - max_cols + 1];
                        // printf("\nZ_row[%ld]=%ld\n", counter_col + Z_row_length - max_cols + 1, Z_row[counter_col + Z_row_length - max_cols + 1]);
                     }
                     else
                     {
                        diag = (isUnknownBase(ctx.X[col]) ? SUBSTITUTION_UNKNOWN_COST : (isSameBase(ctx.X[col], ctx.Y[row]) ? 0 : SUBSTITUTION_COST)) + prev_z;
                        // printf("\nprev_z=%ld\n", prev_z);
                     }

                     long top = INSERTION_COST + Z_row[counter_col + Z_row_length - max_cols];
                     long right = INSERTION_COST + Y_col[row];
                     // printf("Z ROW diag=%ld, top=%ld, right=%ld\n", diag, top, right);
                     // print_array(Z_row, Z_row_length + 1);
                     prev_value = Y_col[row];
                     // printf("ctx.X[%ld]=%c, ctx.Y[%ld]=%c\n", col, ctx.X[col], row, ctx.Y[row]);
                     // printf("Y_col[row]=%ld\n", prev_value);
                     // printf("counter_col=%ld, counter_row=%ld, z_row_val=%ld\n", counter_col, counter_row, Z_row[counter_col]);
                     // print_array(Z_row, Z_row_length + 1);
                     Y_col[row] = min(diag, right, top);
                     // print_array(Y_col, N + 1);
                     // printf("\n");
                  }
                  // print_array(Z_row, Z_row_length + 1);
               }
               value_to_keep = Y_col[row];
            }
            // print_array(Y_col, N + 1);
            // printf("value_to_keep = %ld, prev=%ld\n", value_to_keep, prev_value);
            prev_z = Z_row[counter_col + Z_row_length - max_cols];
            // printf("prev_z=%ld\n", prev_z);
            Z_row[counter_col + Z_row_length - max_cols] = value_to_keep;
            if (counter_col == max_cols - 1)
            {
               // printf("\nput prev_value=%ld in Z_row[%ld]\n\n", prev_value, counter_col + 1);
               Z_row[Z_row_length] = prev_value;
            }
            // printf("put value_to_keep=%ld in Z_row[%ld]\n", value_to_keep, counter_col);
            // printf("\n");
            // print_array(Y_col, N + 1);
            // printf("z_row=%d, z_col=%d\n", z_row, z_col);
            // print_array(Z_row, Z_row_length + 1);
         }
      }
   }
   free(Y_col);
   free(Z_row);
   return Y_col[0];
}
long EditDistance_NW_Iter_A(char *A, size_t lengthA, char *B, size_t lengthB)
{
   _init_base_match();
   struct NW_MemoContext ctx;

   // printf("A = %s", A);
   // printf("B = %s", B);

   if (lengthA >= lengthB)
   {
      ctx.X = A;
      ctx.M = lengthA;
      ctx.Y = B;
      ctx.N = lengthB;
   }
   else
   {
      ctx.X = B;
      ctx.M = lengthB;
      ctx.Y = A;
      ctx.N = lengthA;
   }
   size_t M = ctx.M;
   size_t N = ctx.N;
   long *colonne = malloc((ctx.N + 1) * sizeof(long));

   // Initialisation de la dernière colonne
   colonne[ctx.N] = 0;
   for (int i = ctx.N - 1; i >= 0; i--)
   {
      colonne[i] = 2 * isBase(ctx.Y[i]) + colonne[i + 1];
   }
   long top, right, diag;
   // Remplissage du tableau
   for (int x = ctx.M - 1; x >= 0; x--)
   {
      for (int y = ctx.N; y >= 0; y--)
      {

         if (y == ctx.N)
         {
            diag = colonne[y];
            colonne[y] = (isBase(ctx.X[x]) ? INSERTION_COST : 0) + colonne[y];
         }
         else
         {

            if (!isBase(ctx.X[x]))
            {
               diag = colonne[y];
               ManageBaseError(ctx.X[x]);
            }
            else if (!isBase(ctx.Y[y]))
            {
               diag = colonne[y];
               colonne[y] = colonne[y + 1];
               ManageBaseError(ctx.Y[y]);
            }
            else
            {
               long val_diag = (isUnknownBase(ctx.X[x]) ? SUBSTITUTION_UNKNOWN_COST : (isSameBase(ctx.X[x], ctx.Y[y]) ? 0 : SUBSTITUTION_COST)) + diag;
               long left = INSERTION_COST + colonne[y + 1];
               long top = INSERTION_COST + colonne[y];
               diag = colonne[y];
               colonne[y] = min(val_diag, left, top);
            }
         }
      }
      // On affiche la colonne
      for (int i = 0; i < ctx.N + 1; i++)
      {
         printf("%.2ld ", colonne[i]);
      }
      printf("\n");
   }
   return colonne[0];
}

void calcul_bloc(long *Y_col, long *X_little_row, int deb_x, int deb_y, int K_x, int K_y, struct NW_MemoContext ctx)
{
   long diag = 0;

   // Y_col et X_little_row pas encore calculé
   if ((Y_col[deb_y] == -1) && (X_little_row[deb_x] == -1))
   {

      Y_col[deb_y] = 0;
      diag = 0;
      X_little_row[deb_x] = 0;

      for (int y = deb_y - 1; (y >= deb_y - K_y) && (y >= 0); y--)
      {
         Y_col[y] = (isBase(ctx.Y[y]) ? INSERTION_COST : 0) + Y_col[y + 1];
      }

      for (int x = deb_x - 1; (x >= deb_x - K_x) && (x >= 0); x--)
      {
         X_little_row[x] = (isBase(ctx.X[x]) ? INSERTION_COST : 0) + X_little_row[x + 1];
      }


   }


   // Y_col pas encore calculé
   if (Y_col[deb_y] == -1)
   {
      Y_col[deb_y] = X_little_row[deb_x] + (isBase(ctx.Y[deb_y]) ? INSERTION_COST : 0);
      diag = Y_col[deb_y];
      for (int y = deb_y - 1; (y >= deb_y - K_y) && (y >= 0); y--)
      {
         Y_col[y] = (isBase(ctx.Y[y]) ? INSERTION_COST : 0) + Y_col[y + 1];
      }
   }
   // X_little_row pas encore calculé
   if (X_little_row[deb_x] == -1)
   {
      X_little_row[deb_x] = Y_col[deb_y] + (isBase(ctx.X[deb_x]) ? INSERTION_COST : 0);
      diag = X_little_row[deb_x];
      for (int x = deb_x - 1; (x >= deb_x - K_x) && (x >= 0); x--)
      {
         X_little_row[x] = (isBase(ctx.X[x]) ? INSERTION_COST : 0) + X_little_row[x + 1];
      }
   }
   // Y_col et X_little_row calculé
   for (int x = deb_x; (x >= deb_x - K_x) && (x >= 0); x--)
   {
      for (int y = deb_y; (y >= deb_y - K_y) && (y >= 0); y--)
      {

         if (!isBase(ctx.X[x]))
         {
            diag = Y_col[y];
            ManageBaseError(ctx.X[x]);
         }
         else if (!isBase(ctx.Y[y]))
         {
            diag = Y_col[y];
            Y_col[y] = Y_col[y + 1];
            ManageBaseError(ctx.Y[y]);
         }
         else
         {

            long val_diag = (isUnknownBase(ctx.X[x]) ? SUBSTITUTION_UNKNOWN_COST : (isSameBase(ctx.X[x], ctx.Y[y]) ? 0 : SUBSTITUTION_COST)) + diag;
            long right = INSERTION_COST + Y_col[y];
            long top = INSERTION_COST + Y_col[y + 1];
            diag = Y_col[y];
            Y_col[y] = min(val_diag, right, top);
         }
         X_little_row[x] = Y_col[y];
      }
      // On affiche la colonne
      // for (int i = 0; i < ctx.N + 1; i++)
      // {
      //    printf("%.2ld ", Y_col[i]);
      // }
      // printf("\n");
   }
}

long test_calcul_bloc(char *A, size_t lengthA, char *B, size_t lengthB)
{
   _init_base_match();
   struct NW_MemoContext ctx;

   // printf("A = %s", A);
   // printf("B = %s", B);

   if (lengthA >= lengthB)
   {
      ctx.X = A;
      ctx.M = lengthA;
      ctx.Y = B;
      ctx.N = lengthB;
   }
   else
   {
      ctx.X = B;
      ctx.M = lengthB;
      ctx.Y = A;
      ctx.N = lengthA;
   }
   size_t M = ctx.M;
   size_t N = ctx.N;
   long *colonne = malloc((ctx.N + 1) * sizeof(long));
   long *ligne = malloc((ctx.M + 1) * sizeof(long));
   for (int i = 0; i < (ctx.N + 1); i++)
   {
      colonne[i] = -1;
   }
   for (int i = 0; i < (ctx.M + 1); i++)
   {
      ligne[i] = -1;
   }

   calcul_bloc(colonne, ligne, ctx.M, ctx.N, ctx.M + 1, ctx.N + 1, ctx);

   return colonne[0];
}