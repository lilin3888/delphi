#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>

using namespace std;



struct position
{
    float x,y,z;
};

struct atomInfo
{
    position xyz;
    float charge,radius;
    string atomHead,atomName,ResidueName,chain;
    int atomNumber,residueNumber;
};

class pdbInfo
{
public:
    int NAtom;
    position center;
    atomInfo* atoms;
    void readpdb( string filename)
    {

        this->center.x=1;
        cout << "read pdb: center.x: "<< this->center.x << endl;


        ifstream filein;
        string line;
        int iatom;
        filein.open( filename.c_str() );

        //###################### get the atom onumber: #################
        this->NAtom=0;
        while(!filein.eof())
        {

            getline(filein,line);

            //cout << "[ " << line << " ]" << endl;

            if(line.substr(0,6)=="ATOM  " || line.substr(0,6)=="HETATM")
            {
                //cout << line.substr(0,6) << endl;
                //cout << line << endl;

                this->NAtom++;
            }

        }
        cout << "Atom number is: " << this->NAtom << endl;

        filein.close();



        atomInfo * atomA = new atomInfo [this->NAtom];

        //###################### Read file: #################
        iatom=-1;


        filein.open("1ARB.pqr");
        while(!filein.eof())
        {

            getline(filein,line);

            //cout << "[ " << line << " ]" << endl;

            if(line.substr(0,6)=="ATOM  " || line.substr(0,6)=="HETATM")
            {
                iatom++;
                atomA[iatom].atomHead = line.substr(0,6);
                atomA[iatom].atomNumber = atoi(line.substr(6,5).c_str());
                atomA[iatom].atomName = line.substr(12,4);
                atomA[iatom].ResidueName = line.substr(17,3);
                atomA[iatom].chain = line.substr(21,1);

                atomA[iatom].xyz.x = atof(line.substr(30,8).c_str());
                atomA[iatom].xyz.y = atof(line.substr(38,8).c_str());
                atomA[iatom].xyz.z = atof(line.substr(46,8).c_str());
                atomA[iatom].charge = atof(line.substr(54,8).c_str());
                atomA[iatom].radius = atof(line.substr(62,7).c_str());


                cout << showpoint << fixed;
                cout << atomA[iatom].atomHead
                     << setw(5)<< atomA[iatom].atomNumber
                     << " " << atomA[iatom].atomName
                     << " " << atomA[iatom].ResidueName
                     << " " << atomA[iatom].chain
                     << setprecision(3) << setw(8) << atomA[iatom].xyz.x
                     << setprecision(3) << setw(8) << atomA[iatom].xyz.y
                     << setprecision(3) << setw(8) << atomA[iatom].xyz.z
                     << setprecision(4) << setw(8) << atomA[iatom].charge
                     << setprecision(4) << setw(7) << atomA[iatom].radius
                     << endl;
            }


        }

        filein.close();

    }


};

int main(int argc, char *argv[])
{
    int i,iatom;
    int iAtomNum;
    string line;
/*
    cout << "argc: " << argc << endl;
    if(argc > 1)
    {
        cout << "argv[1]: " << argv[1] << endl;
    }

    ifstream filein;
    filein.open("1ARB.pqr");

    //###################### get the atom onumber: #################
    iAtomNum=0;
    while(!filein.eof())
    {

        getline(filein,line);

        //cout << "[ " << line << " ]" << endl;

        if(line.substr(0,6)=="ATOM  " || line.substr(0,6)=="HETATM")
        {
            //cout << line.substr(0,6) << endl;
            //cout << line << endl;

            iAtomNum++;
        }

    }
    cout << "Atom number is: " << iAtomNum << endl;

    filein.close();



    atomInfo * atomA = new atomInfo [iAtomNum];

    //###################### Read file: #################
    iatom=-1;


    filein.open("1ARB.pqr");
    while(!filein.eof())
    {

        getline(filein,line);

        //cout << "[ " << line << " ]" << endl;

        if(line.substr(0,6)=="ATOM  " || line.substr(0,6)=="HETATM")
        {
            iatom++;
            atomA[iatom].atomHead = line.substr(0,6);
            atomA[iatom].atomNumber = atoi(line.substr(6,5).c_str());
            atomA[iatom].atomName = line.substr(12,4);
            atomA[iatom].ResidueName = line.substr(17,3);
            atomA[iatom].chain = line.substr(21,1);

            atomA[iatom].xyz.x = atof(line.substr(30,8).c_str());
            atomA[iatom].xyz.y = atof(line.substr(38,8).c_str());
            atomA[iatom].xyz.z = atof(line.substr(46,8).c_str());
            atomA[iatom].charge = atof(line.substr(54,8).c_str());
            atomA[iatom].radius = atof(line.substr(62,7).c_str());

            cout << showpoint << fixed;
            cout << atomA[iatom].atomHead
                 << setw(5)<< atomA[iatom].atomNumber
                 << " " << atomA[iatom].atomName
                 << " " << atomA[iatom].ResidueName
                 << " " << atomA[iatom].chain
                 << setprecision(3) << setw(8) << atomA[iatom].xyz.x
                 << setprecision(3) << setw(8) << atomA[iatom].xyz.y
                 << setprecision(3) << setw(8) << atomA[iatom].xyz.z
                 << setprecision(4) << setw(8) << atomA[iatom].charge
                 << setprecision(4) << setw(7) << atomA[iatom].radius
                 << endl;
        }
    }

    filein.close();

    */


    //####copy pdb#####################################################
/*
    atomInfo * atomB = new atomInfo [iAtomNum];

    for (i=0; i<iAtomNum; i++)
    {
        atomB[i]=atomA[i];
    }

    pdbInfo pdbA;
    pdbA.readpdb("1ARB.pqr");
    cout << "pdbA.center: " << pdbA.center.x << endl;

    pdbInfo pdbB;
    pdbB=pdbA;
    cout << "pdbB.center: " << pdbB.center.x << endl;


    delete [] atomA;

*/
    pdbInfo pdbA;

    pdbA.readpdb("1ARB.pqr");




}
