#! /bin/csh

foreach i (`ls *f90`)
echo '--------------------------'
echo $i
diff $i ../source/$i

end
