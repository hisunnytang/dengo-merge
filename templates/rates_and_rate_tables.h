/*

The generalized rate data type holders.

*/

typedef struct RateTable {
  int nbins;
  /* Note that we could include a piece of data that points to the 
     variable on which this rate depends; but we do not. */
  double dbin;
  double idbin;
  double bounds[2];
  double *values;
} RateTable;

typedef struct TableOfRates {
  int nrates;
  RateTable *rates;
} TableOfRates;
