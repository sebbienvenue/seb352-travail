Readme for the source of RuNNer
-------------------------------

Programming guidelines:

1) Use unique variable names!
   The only variable names with multiple use (and maybe dimension) are
   - i1,i2,i3,i4 etc
   - arrays ending with _local, _mpi or _temp

2) Reuse subroutines for multiple purposes if possible to avoid redundant source code
   Avoid array names which are used in a non-standard way (use _local, _mpi or _temp instead inside)

3) comment what you do, in particular if you are not sure if it is right or working

4) In a call statement sort arguments: integers, real*8, logicals, characters

5) When passing data in a call statement to subroutines, keep the following order:
First all integers, then real*8, then character, finally logicals.

6) use insets for if, do etc.

7) Everything that is not yet working or might be wrong in the code must be labeled by 'FIXME' and a comment text 

-----------------------------------------------------------------------------------------------

Todo:
- think again about Kalman parameters, are they useful for large data sets?

- why does fitting with worst 10 % of energies have strongly varying numbers of points in each epoch???

- would it be useful to distribute only the initial bias weights in a special way? Bias weights should be at results of linear combinations. Probably at the moment they are very far and at too low values.

- change definitions of angles in angular pair symmetry functions

- improve preconditioner

- updatebyelement requires revision for the forces: At the moment, if only weights of element A are updated, then the force RMSE is calculated only for atoms of this element, but also the forces on atoms of other elements depend in principle on the weights of element A, (because A can be in their neighborhood). Further, even for the force of an atom of element A only the weights of this element A are updated, but the force on A depends also on the weights for the other elements, because they can be in the environment. So this is all not logic.


