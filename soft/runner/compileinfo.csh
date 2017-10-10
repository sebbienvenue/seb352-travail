#! /bin/csh

rm -f compileinfomiddle
touch compileinfomiddle
set i = `hostname`
set j = `whoami`
set k = `date`
set r = `svn info | grep Revision | cut -c 11-16`
echo '      write(ounit,*)"Compiled by user ' $j ' on host ' $i '"' >> compileinfomiddle
echo '      write(ounit,*)"Compilation time ' $k '"' >> compileinfomiddle
echo '      write(ounit,*)"Source Code revision ' $r '"' >> compileinfomiddle

