for pn in $(cat list.txt)
do 
	pn2=$(echo $pn|sed 's/f95/cpp/g;s/inc_epsmakmod_//g')
	echo $pn $pn2
	
	./f2c.sh ./space_module_fortran2/$pn ./auto_cpp/$pn2

done
