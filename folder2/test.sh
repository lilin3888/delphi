while read line
do
    flag=0
    string=$line


    flag=$(echo $string | grep "select case" )   #>> $1.cpp

#    echo "flag:",$flag

    if [[ $flag != "" ]]
    then
        #echo $string
        echo $string|sed 's/select case/swich /g;s/;/{/g' 
    elif [[ ${string:0:4} == "case" ]]
    then
	echo $string | sed 's/(/ /;s/)//'
    else
        pn=1
	#echo $string
    fi
done < f2c_v.tmp 

