if [[ $#  != 2 ]]
then
	echo "there are totally " $# "parameters:" $1 
	echo "stop: 2 parameters are expected."
	exit 0
fi

###### replace variables: #######
cp $1 f2c_v.tmp
echo $1
#head $1
while read line
do
	string=$line
	flag=${string:0:1}
	if [[ $flag == "!" ]] 
	then
		continue
	fi

	pn1=$(echo $string |awk '{print $1}')
	pn2=$(echo $string |awk '{print $2}')
#	echo $pn1 $pn2
	
#	echo "before grep..."
	sed "s/$pn1/$pn2/g" f2c_v.tmp > f2c_v2.tmp
#	echo "after grep"
	
#	echo ""
#	head f2c_v.tmp

	cp f2c_v2.tmp f2c_v.tmp

done < variables.list 

###### replace operators: ######
sed     's/type(coord) ::/SGrid <float>/' f2c_v.tmp | \
sed     's/type(int_coord) ::/SGrid <int>/ ' | \
sed     's/implicit none//' |\
sed     's/subroutine/void/' |\
sed     's/logical ::/bool/' |\
sed     's/integer ::/int/' |\
sed     's/real ::/float/' |\
sed     's/logical /bool /' |\
sed     's/integer /int /' |\
sed     's/real /float /' |\
sed     's/real(/Int2Float(/' |\
sed     's/float(/Int2Float(/' |\
sed     's/int(/Float2Int(/' |\
sed     's/\!/\/\//' |\
sed     's/\;/\n/g' |\
#sed     's/\&/\\/' |\
sed     's/\.true\./true/g' |\
sed     's/\.false\./false/g' |\
sed     's/\.eq\./==/g' |\
sed     's/\.ne\./!=/g' |\
sed     's/\.lt\./</g' |\
sed     's/\.gt\./>/g' |\
sed     's/\.le\./<=/g' |\
sed     's/\.ge\./>=/g' |\
sed     's/\.not\./\!/g' |\
sed     's/\.and\./\&\&/g' |\
sed     's/\.or\./\|\|/g'|\
#sed    's/do /for /' |\
sed     's/end/\}\/\// ' |\
sed     's/then/\{/' > f2c_v2.tmp
#sed    's///' |\
#sed    's///' |\

cp f2c_v2.tmp f2c_v.tmp

###### some complicated operations:######

###### do for ######
while read line
do
    flag=0
    string=$line


    flag=$(echo $string|awk '{if($1 == "do" ) print "do"}')   #>> $1.cpp


#    echo "flag:",$flag

    if [[ $flag == "do" ]]
    then
#       echo $flag $string
        index=$(echo $string | sed 's/=/ /;s/,/ /'|awk '{print $2}' )
        start=$(echo $string | sed 's/=/ /;s/,/ /'|awk '{print $3}' )
        finish=$(echo $string | sed 's/=/ /;s/,/ /'|awk '{print $4}' )
        #echo $string
        #echo $index $start $finish
        echo $index $start $finish|awk '{printf("for(%s=%s;%s<=%s;%s++){\n",$1,$2,$1,$3,$1)}'
#    elif
#       echo "elif"
    else 
        #echo "else"
        echo $string
    fi

done < f2c_v.tmp > f2c_v2.tmp

cp f2c_v2.tmp f2c_v.tmp

###### add ; ######

awk '{if(substr($1,1,2) != "//" && $0 != "") {printf("%s%s\n",$0,";")} 
else print $0}' f2c_v.tmp > f2c_v2.tmp

cp f2c_v2.tmp f2c_v.tmp

sed 's/{;/{/g; s/else;/}\nelse/g; s/else/else{/g;s/if;/if/g; s/do;/do/g' f2c_v.tmp > f2c_v2.tmp

cp f2c_v2.tmp f2c_v.tmp

###### remove & & continue line ######

awk '{if(substr($0,length($0)-1,2)=="&;") printf("%s",substr($0,1,length($0)-2)) 
	else if(substr($0,1,1)=="&") print substr($0,2,length($0)-1)
	else print $0}' f2c_v.tmp > f2c_v2.tmp


cp f2c_v2.tmp f2c_v.tmp

###### write and print ######
while read line
do
    flag=0
    string=$line


    flag=$(echo $string | grep "write(6,\*)\|print" )   #>> $1.cpp


#    echo "flag:",$flag

    if [[ $flag != "" ]]
    then
        #echo $string
        echo $string|sed 's/write(6,\*)/cout <</g;s/,/ << /g;s/;/ << endl;/g' | sed s/\'/\"/g |sed 's/print \*/cout <</'
    else
        echo $string
    fi
done < f2c_v.tmp > f2c_v2.tmp

cp f2c_v2.tmp f2c_v.tmp

###### case select#####
while read line
do
    flag=0
    string=$line

    flag=$(echo $string | grep "select case" )   #>> $1.cpp

#    echo "flag:",$flag

    if [[ $flag != "" ]]
    then
        #echo $string
        echo $string|sed 's/select case/switch /g;s/;/{/g'
    elif [[ ${string:0:4} == "case" ]]
    then
        echo $string | sed 's/(/ /;s/)//;s/;/:/;s/case default/default/'
    else
        echo $string
    fi
done < f2c_v.tmp > f2c_v2.tmp


cp f2c_v2.tmp f2c_v.tmp


###### handle arrays ######

for pn in $(grep "^int \|^bool \|^float \|^SGrid " f2c_v2.tmp |sed 's/,/ /g')
do 
        echo $pn 
done|grep "("|sed 's/(/ /g'|awk '{print $1}' > array.lst

while read line
do
    string=$line 
    flag_g=0
    for pn in $(cat array*.lst)
    do
        flag=$(echo $string |grep "$pn(" )

        if [[ $flag != ""  ]]
        then

        flag_g=1
        pn2=${pn}[
        len=${#pn2}
        #echo $len
        newline=$(echo $string|sed "s/$pn(/$pn[/g"|awk -v str=$pn2 -v len=$len '{L=length($0);flag=-100;flag2=-100}{for(i=1;i<=L;i++){a=substr($0,i,1);b=substr($0,i,len);if(a=="(") flag=flag+1; if(b==str) {flag2=1; flag=1};if(a==")") flag=flag-1; if(a==")"&& flag==0 && flag2==1){printf("]") ;flag2=0} else printf("%s",a)} }END{printf("\n")}')
#       echo $newline
        string=$newline
        fi

    done

    flag_b=$(echo $string|awk '{print $1}')
    if [[ $flag_g == "1" && $flag_b != "bool" && $flag_b != "int" && $flag_b != "float" && $flag_b != "SGrid" ]]
    then
        echo $string | sed 's/,/][/g'
    else
        echo $string
    fi
done< f2c_v.tmp > f2c_v2.tmp


cp f2c_v2.tmp f2c_v.tmp


###### add header ######


cat header.txt f2c_v.tmp > f2c_v2.tmp


cp f2c_v2.tmp f2c_v.tmp


cp f2c_v.tmp $2
