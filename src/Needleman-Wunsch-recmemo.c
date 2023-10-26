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
      printf("%ld, ", array[i]);
   }
   printf("]\n");
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

   long Y_col_size = N + 1 > 4096 ? 4096 : N + 1;
   // long Y_col_size = 4096;

   long *Y_col = (long *)malloc(Y_col_size * sizeof(long));
   if (Y_col == NULL)
   {
      perror("EditDistance_NW_Iter: malloc of Y_col. You must have at least one!");
      exit(EXIT_FAILURE);
   }
   Y_col[N] = 0;
   for (int row = N + 1; row >= N + 1 - Y_col_size; row--)
   {
      Y_col[row] = (isBase(ctx.Y[row]) ? INSERTION_COST : 0) + Y_col[row + 1];
   }

   long prev_value;

   for (int col = M - 1; col >= 0; col--)
   {
      for (int row = N; row >= 0; row--)
      {
         if (row % Y_col_size == 0)
         {
         }
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
            long left = INSERTION_COST + Y_col[row + 1];
            long top = INSERTION_COST + Y_col[row];
            prev_value = Y_col[row];
            Y_col[row] = min(diag, left, top);
         }
      }
   }

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
   long *colonne = malloc(ctx.N * sizeof(long));
   long *ligne = malloc(ctx.M * sizeof(long));
   // Initialisation de la dernière ligne
   ligne[ctx.M - 1] = 0;
   for (int i = ctx.M - 2; i >= 0; i--)
   {
      ligne[i] = 2 * isBase(ctx.X[i]) + ligne[i + 1];
   }
   // Initialisation de la dernière colonne
   colonne[ctx.N - 1] = 0;
   for (int i = ctx.N - 2; i >= 0; i--)
   {
      colonne[i] = 2 * isBase(ctx.Y[i]) + colonne[i + 1];
   }
   // Remplissage du tableau
   for (int i = ctx.N - 2; i >= 0; --i)
   {
      for (int j = ctx.M - 2; j >= 0; --j)
      {
         if (!isBase(ctx.X[j]))
         {
            ligne[j] = ligne[j + 1];
         }
         if (!isBase(ctx.Y[i]))
         {
            ligne[j] = ligne[j];
         }
         if (isBase(ctx.X[j]) && isBase(ctx.Y[i]))
         {
            long diag = (isUnknownBase(ctx.X[j]) ? SUBSTITUTION_UNKNOWN_COST : (isSameBase(ctx.X[j], ctx.Y[i]) ? 0 : SUBSTITUTION_COST)) + ligne[j + 1];
            long left = INSERTION_COST + ligne[j];
            long top = INSERTION_COST + colonne[i];
            ligne[j] = min(diag, left, top);
         }
      }
   }
}
