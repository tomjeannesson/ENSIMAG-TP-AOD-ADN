/**
 * \file Needleman-Wunsch-recmemo.c
 * \brief recursive implementation with memoikation of Needleman-Wunsch global alignment algorithm that computes the distance between two genetic sequences
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

/* Context of the memoikation : passed to all recursive calls */
/** \def NOT_YET_COMPUTED
 * \brief default value for memoikation of minimal distance (defined as an impossible value for a distance, -1).
 */
#define NOT_YET_COMPUTED -1L

/** \struct NW_MemoContext
 * \brief data for memoikation of recursive Needleman-Wunsch algorithm
 */
struct NW_MemoContext
{
   char *X;     /*!< the longest genetic sequences */
   char *Y;     /*!< the shortest genetic sequences */
   size_t M;    /*!< length of X */
   size_t N;    /*!< length of Y,  N <= M */
   long **memo; /*!< memoikation table to store memo[0..M][0..N] (including stopping conditions phi(M,j) and phi(i,N) */
};

/*
 *  static long EditDistance_NW_RecMemo(struct NW_MemoContext *c, size_t i, size_t j)
 * \brief  EditDistance_NW_RecMemo :  Private (static)  recursive function with memoikation \
 * direct implementation of Needleman-Wursch extended to manage FASTA sequences (cf TP description)
 * \param c : data passed for recursive calls that includes the memoikation array
 * \param i : starting position of the left sequence :  c->X[ i .. c->M ]
 * \param j : starting position of the right sequence :  c->Y[ j .. c->N ]
 */
static long EditDistance_NW_RecMemo(struct NW_MemoContext *c, size_t i, size_t j)
/* compute and returns phi(i,j) using data in c -allocated and initialiked by EditDistance_NW_Rec */
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
         long min = /* initialikation  with cas 1*/
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
 * It allocates and initailikes data (NW_MemoContext) for memoikation and call the
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
   { /* Allocation and initialikation of ctx.memo to NOT_YET_COMPUTED*/
      /* Note: memo is of sike (N+1)*(M+1) but is stored as (M+1) distinct arrays each with (N+1) continuous elements
       * It would have been possible to allocate only one big array memekone of (M+1)*(N+1) elements
       * and then memo as an array of (M+1) pointers, the memo[i] being the address of memkone[i*(N+1)].
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
      }
   }
   long res = Y_col[0];
   free(Y_col);
   return res;
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
   long K = 100;
   long number_of_K_rows = N / K + 1;
   long number_of_K_cols = M / K + 1;
   long K_row_length = K > M ? M : K;
   long *K_row = (long *)malloc((K_row_length + 1) * sizeof(long));
   long prev_k;
   long diag;
   long right;
   long top;
   long max_rows;
   long max_cols;
   long value_to_keep;
   long col;
   long row;
   for (int k_col = number_of_K_cols - 1; k_col >= 0; k_col--)
   {
      max_cols = M - (number_of_K_cols - 1 - k_col) * K > K ? K : M - (number_of_K_cols - 1 - k_col) * K;
      for (int k_row = number_of_K_rows - 1; k_row >= 0; k_row--)
      {
         max_rows = N - (number_of_K_rows - 1 - k_row) * K > K ? K : N - (number_of_K_rows - 1 - k_row) * K;
         if (k_row != number_of_K_rows - 1)
         {
            max_rows--;
         }
         for (int counter_col = max_cols - 1; counter_col >= 0; counter_col--)
         {
            col = k_col > 0 ? counter_col + (k_col - 1) * K + M - (number_of_K_cols - 1) * K : counter_col;
            for (int counter_row = max_rows; counter_row >= 0; counter_row--)
            {
               row = k_row > 0 ? counter_row + (k_row - 1) * K + N - (number_of_K_rows - 1) * K : counter_row;

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
                     Y_col[row] = K_row[counter_col + K_row_length - max_cols];
                  }
                  else
                  {

                     Y_col[row] = Y_col[row + 1];
                  }
                  ManageBaseError(ctx.Y[row]);
               }
               else
               {
                  if (k_row == number_of_K_rows - 1 || counter_row != max_rows)
                  {

                     diag = (isUnknownBase(ctx.X[col]) ? SUBSTITUTION_UNKNOWN_COST : (isSameBase(ctx.X[col], ctx.Y[row]) ? 0 : SUBSTITUTION_COST)) + prev_value;
                     top = INSERTION_COST + Y_col[row + 1];
                     right = INSERTION_COST + Y_col[row];
                     prev_value = Y_col[row];
                     Y_col[row] = min(diag, right, top);
                  }
                  else
                  {
                     diag;
                     if (counter_col == max_cols - 1)
                     {

                        diag = (isUnknownBase(ctx.X[col]) ? SUBSTITUTION_UNKNOWN_COST : (isSameBase(ctx.X[col], ctx.Y[row]) ? 0 : SUBSTITUTION_COST)) + K_row[counter_col + K_row_length - max_cols + 1];
                     }
                     else
                     {
                        diag = (isUnknownBase(ctx.X[col]) ? SUBSTITUTION_UNKNOWN_COST : (isSameBase(ctx.X[col], ctx.Y[row]) ? 0 : SUBSTITUTION_COST)) + prev_k;
                     }

                     top = INSERTION_COST + K_row[counter_col + K_row_length - max_cols];
                     right = INSERTION_COST + Y_col[row];
                     prev_value = Y_col[row];
                     Y_col[row] = min(diag, right, top);
                  }
               }
               value_to_keep = Y_col[row];
            }
            prev_k = K_row[counter_col + K_row_length - max_cols];
            K_row[counter_col + K_row_length - max_cols] = value_to_keep;
            if (counter_col == max_cols - 1)
            {
               K_row[K_row_length] = prev_value;
            }
         }
      }
   }
   long res = Y_col[0];
   free(Y_col);
   free(K_row);
   return res;
}

long CO_BREAKPOINT = 10;

long CalculateBlock_CO(struct NW_MemoContext ctx, long *X_row, long X_start, long X_end, long *Y_col, long Y_start, long Y_end)
{
   long prev_value_y;
   long prev_value_k;
   long value_to_keep;
   long diag;
   long left;
   long top;

   if (X_end - X_start > CO_BREAKPOINT)
   {
      CalculateBlock_CO(ctx, X_row, X_start + (X_end - X_start) / 2 + 1, X_end, Y_col, Y_start, Y_end);
      CalculateBlock_CO(ctx, X_row, X_start, X_start + (X_end - X_start) / 2, Y_col, Y_start, Y_end);
   }
   else if (Y_end - Y_start > CO_BREAKPOINT)
   {
      CalculateBlock_CO(ctx, X_row, X_start, X_end, Y_col, Y_start + (Y_end - Y_start) / 2 + 1, Y_end);
      CalculateBlock_CO(ctx, X_row, X_start, X_end, Y_col, Y_start, Y_start + (Y_end - Y_start) / 2);
   }
   else
   {

      for (int col = X_end; col >= X_start; col--)
      {
         for (int row = Y_end; row >= Y_start; row--)
         {
            if (row == ctx.N)
            {
               prev_value_y = Y_col[row];
               Y_col[row] = (isBase(ctx.X[col]) ? INSERTION_COST : 0) + Y_col[row];
            }
            else if (!isBase(ctx.X[col]))
            {
               prev_value_y = Y_col[row];
               ManageBaseError(ctx.X[col]);
            }
            else if (!isBase(ctx.Y[row]))
            {
               prev_value_y = Y_col[row];
               if (row == Y_end)
               {
                  Y_col[row] = X_row[col];
               }
               else
               {
                  Y_col[row] = Y_col[row + 1];
               }
               ManageBaseError(ctx.Y[row]);
            }
            else
            {
               if (row == Y_end)
               {
                  if (col == X_end)
                  {

                     diag = (isUnknownBase(ctx.X[col]) ? SUBSTITUTION_UNKNOWN_COST : (isSameBase(ctx.X[col], ctx.Y[row]) ? 0 : SUBSTITUTION_COST)) + X_row[col + 1];
                  }
                  else
                  {
                     diag = (isUnknownBase(ctx.X[col]) ? SUBSTITUTION_UNKNOWN_COST : (isSameBase(ctx.X[col], ctx.Y[row]) ? 0 : SUBSTITUTION_COST)) + prev_value_k;
                  }
                  left = INSERTION_COST + Y_col[row];
                  top = INSERTION_COST + X_row[col];
               }
               else
               {
                  diag = (isUnknownBase(ctx.X[col]) ? SUBSTITUTION_UNKNOWN_COST : (isSameBase(ctx.X[col], ctx.Y[row]) ? 0 : SUBSTITUTION_COST)) + prev_value_y;
                  left = INSERTION_COST + Y_col[row];
                  top = INSERTION_COST + Y_col[row + 1];
               }
               prev_value_y = Y_col[row];
               Y_col[row] = min(diag, left, top);
            }
            value_to_keep = Y_col[row];
         }
         prev_value_k = X_row[col];
         X_row[col] = value_to_keep;
         if (col == X_end)
         {
            X_row[X_end + 1] = prev_value_y;
         }
      }
   }
}

long EditDistance_NW_Iter_CO(char *A, size_t lengthA, char *B, size_t lengthB)
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

   long *X_row = (long *)malloc((M + 1) * sizeof(long));
   X_row[M] = 0;
   for (int row = M - 1; row >= 0; row--)
   {
      X_row[row] = (isBase(ctx.X[row]) ? INSERTION_COST : 0) + X_row[row + 1];
   }

   CalculateBlock_CO(ctx, X_row, 0, M - 1, Y_col, 0, N);

   long res = Y_col[0];
   free(Y_col);
   free(X_row);
   return res;
}