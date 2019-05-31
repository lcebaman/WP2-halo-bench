#!/bin/sh

rm total.txt computation.txt communication.txt

for file in *out
do
    
    for i in 2 
    do 
	echo  $file | awk -F'[_]' ' {printf "%s  ", $2}' >>total.txt
	
	awk '{printf "%f\n",$'${i}'}' ${file} | sort -n  |
	
	awk '
  BEGIN {
    c = 0;
    sum = 0;
  }
  $1 ~ /^[0-9]*(\.[0-9]*)?$/ {
    a[c++] = $1;
    sum += $1;
    sumsq[c] =$1;
  }
  END {
    ave = sum / c;
    if( (c % 2) == 1 ) {
      median = a[ int(c/2) ];
    } else {
      median = ( a[c/2] + a[c/2-1] ) / 2;
    }
    for(x=1; x<=c; x++)
      sq += (sumsq[x]-(sum/c))^2

    OFS="\t";
    print ave, median, a[0], a[c-1], sqrt(sq/c);
  }
'  >>total.txt
     sort -n -k1,1 -o total.txt total.txt

    done

    for i in 3
    do 
	
	echo  $file | awk -F'[_]' ' {printf "%s  ", $2}'>>computation.txt
	
	awk '{printf "%f\n",$'${i}'}' ${file} | sort -n  |
	
	awk '
  BEGIN {
    c = 0;
    sum = 0;
  }
  $1 ~ /^[0-9]*(\.[0-9]*)?$/ {
    a[c++] = $1;
    sum += $1;
    sumsq[c] =$1;
  }
  END {
    ave = sum / c;
    if( (c % 2) == 1 ) {
      median = a[ int(c/2) ];
    } else {
      median = ( a[c/2] + a[c/2-1] ) / 2;
    }
    for(x=1; x<=c; x++)
      sq += (sumsq[x]-(sum/c))^2

    OFS="\t";
    print ave, median, a[0], a[c-1], sqrt(sq/c);
  }'   >> computation.txt
	sort -n -k1,1 -o computation.txt computation.txt
    done
    
    
    
    for i in 4
    do 
	echo  $file | awk -F'[_]' ' {printf "%s  ", $2}' >>communication.txt
	
	awk '{printf "%f\n",$'${i}'}' ${file} | sort -n  |
	
	awk '
  BEGIN {
    c = 0;
    sum = 0;
  }
  $1 ~ /^[0-9]*(\.[0-9]*)?$/ {
    a[c++] = $1;
    sum += $1;
    sumsq[c] =$1;
  }
  END {
    ave = sum / c;
    if( (c % 2) == 1 ) {
      median = a[ int(c/2) ];
    } else {
      median = ( a[c/2] + a[c/2-1] ) / 2;
    }
    for(x=1; x<=c; x++)
      sq += (sumsq[x]-(sum/c))^2

    OFS="\t";
    print ave, median, a[0], a[c-1], sqrt(sq/c);
  }' >> communication.txt
	sort -n -k 1,1 -o  communication.txt communication.txt
    done
done
