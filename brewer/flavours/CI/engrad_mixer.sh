#!/bin/bash

basis=$1
functional=$2
iroot=$3
number=$4

if [ $iroot -ne 0 ]
then
    cat > engrad$number << !
! $basis $functional engrad
    
%tddft
 nroots 10
 iroot $iroot
 tda TRUE
end
    
* xyzfile 0 1 geom.xyz 
!
else
    cat > engrad$number << !
! $basis $functional engrad

* xyzfile 0 1 geom.xyz
!
fi

