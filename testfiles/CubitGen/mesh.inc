********************************** N O D E S **********************************
*NODE, NSET=ALLNODES
       1,    2.000000e+00,   -1.000000e+00,    2.500000e+00
       2,    2.000000e+00,    1.000000e+00,    2.500000e+00
       3,   -2.000000e+00,    1.000000e+00,    2.500000e+00
       4,   -2.000000e+00,   -1.000000e+00,    2.500000e+00
       5,    2.000000e+00,    1.000000e+00,   -2.500000e+00
       6,    2.000000e+00,   -1.000000e+00,   -2.500000e+00
       7,   -2.000000e+00,   -1.000000e+00,   -2.500000e+00
       8,   -2.000000e+00,    1.000000e+00,   -2.500000e+00
       9,    2.000000e+00,    3.000000e+00,    2.500000e+00
      10,   -2.000000e+00,    3.000000e+00,    2.500000e+00
      11,    2.000000e+00,    3.000000e+00,   -2.500000e+00
      12,   -2.000000e+00,    3.000000e+00,   -2.500000e+00
      13,    2.000000e+00,    5.000000e+00,    2.500000e+00
      14,   -2.000000e+00,    5.000000e+00,    2.500000e+00
      15,    2.000000e+00,    5.000000e+00,   -2.500000e+00
      16,   -2.000000e+00,    5.000000e+00,   -2.500000e+00
      17,    2.000000e+00,    7.000000e+00,    2.500000e+00
      18,   -2.000000e+00,    7.000000e+00,    2.500000e+00
      19,    2.000000e+00,    7.000000e+00,   -2.500000e+00
      20,   -2.000000e+00,    7.000000e+00,   -2.500000e+00
********************************** E L E M E N T S ****************************
*ELEMENT, TYPE=C3D8R, ELSET=block1
       1,       4,       7,       8,       3,       1,       6,       5,       2
       2,       5,      11,      12,       8,       2,       9,      10,       3
*ELEMENT, TYPE=C3D8R, ELSET=block2
       3,      11,      15,      16,      12,       9,      13,      14,      10
       4,      15,      19,      20,      16,      13,      17,      18,      14
********************************** N O D E S E T S **********************************
*NSET, NSET=n_bot
       1,       4,       6,       7,
*NSET, NSET=n_top
      17,      18,      19,      20,