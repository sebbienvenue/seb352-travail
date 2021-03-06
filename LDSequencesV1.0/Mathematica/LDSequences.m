(***********************************************************************
This file was generated automatically by the Mathematica front end.
It contains Initialization cells from a Notebook file, which typically
will have the same name as this file except ending in ".nb" instead of
".m".

This file is intended to be loaded into the Mathematica kernel using
the package loading commands Get or Needs.  Doing so is equivalent to
using the Evaluate Initialization Cells menu command in the front end.

DO NOT EDIT THIS FILE.  This entire file is regenerated automatically 
each time the parent Notebook file is saved in the Mathematica front end.
Any changes you make to this file will be overwritten.
***********************************************************************)







BeginPackage["LDSequences`Halton`"]



\!\(Halton::"\<usage\>" = \*"\"\<Halton[n, base] returns the n-th element of \
the van der Corput sequence in base.\nHalton[n, {\!\(b\_1\), \!\(b\_2\), \
...}] returns the n-th element of the Halton sequence in the bases \
{\!\(b\_1\), \!\(b\_2\), ...}.\>\""\)



Begin["Private`"]



Halton[n_Integer?NonNegative, b_Integer?Positive]:=
  
  FromDigits[ {Reverse[IntegerDigits[n, b] ], 0 }, b]

Halton[n_Integer?NonNegative, b:{___Integer?Positive}]:=
    
    FromDigits[ {Reverse[IntegerDigits[n, #] ], 0 }, #]&/@b;

End[];
EndPackage[]



















BeginPackage["LDSequences`TSNet`", "NumberTheory`NumberTheoryFunctions`"]



TSNet::"usage"=
    "TSNet[n, base, dim] returns the n-th element of the (0,s)-Netz \
construction using hyperderivatives in base.
TSNet[n,dim] is equal to TSNet[n, NextPrime[dim], dim].";

Begin["Private`"]



TSNet[n_Integer?NonNegative, dim_Integer?Positive]:=
    TSNet[n, NextPrime[dim-0.5], dim];

TSNet[n_Integer?NonNegative, base_Integer?Positive, dim_Integer?Positive]:=
    TSNet[n,base, dim, Range[0,dim-1]];

\!\(\(TSNet[n_Integer?NonNegative, \ base_Integer?Positive, \ 
        dim_Integer?Positive, \ bases : {___Integer?NonNegative}] := 
      Module[{coeffs = Reverse@IntegerDigits[n, base], \ CreateCoeffs, j, \ 
          k, cfs}, \[IndentingNewLine]CreateCoeffs = 
          Evaluate /@ \((Table[
                  Table[Binomial[k, \ j - 1]*#\^\(k - j + 1\), \ {k, 0, \ 
                        Length[coeffs] - 1}]\  . \ coeffs, \ {j, \ 1, 
                    Length[coeffs]}] &)\); \[IndentingNewLine]cfs = 
          Mod[\((CreateCoeffs /@ bases)\), \ 
            base]; \[IndentingNewLine]\(FromDigits[{#, 0}, base] &\) /@ 
          cfs\[IndentingNewLine]];\)\)

End[];
EndPackage[]





























































































































































BeginPackage["LDSequences`nAlpha`"];

\!\(\(nAlpha::"\<usage\>" = \*"\"\<nAlpha[n, \[Alpha]] creates the n-th \
element of the one-dimensional {n\[Alpha]} sequence. nAlpha[n, \
{\!\(\[Alpha]\_1\), \!\(\[Alpha]\_2\), ...}] creates the n-th element of the \
multidimensional {n\[Alpha]} sequence using the \!\(\[Alpha]\_i\) given. They \
can e.g. be created using CreateAlphas[n].\>\"";\)\[IndentingNewLine]
  \(CreateAlphas::"\<usage\>" = "\<CreateAlphas[n] creates a list of the Real \
numbers represented by the continued fractions [i,i,i,i,i,i,i,i,...], where 1\
\[LessEqual]i\[LessEqual]n.\>";\)\)

Begin["Private`"]

nAlpha[n_Integer?NonNegative, a_Real?Positive]:=FractionalPart[n*a];
nAlpha[n_Integer?NonNegative, a:{___Real?Positive}]:=nAlpha[n, #]&/@a;

CreateAlphas[n_Integer?Positive]:=
  N[FromContinuedFraction[Table[#, {500}]]&/@Range[n], 50]

End[];
EndPackage[];































