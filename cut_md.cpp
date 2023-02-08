//Cut_md can run classical molecular dynamic using a lennard jones potential to compute forces and energies.
//The velocity Verlet algorithm is employed for solving the Newton's equations of motion.
//INPUT FILE: a list of atoms, one for line, with name and xyz coordinates. No header and no comment
//OUTPUT FILE: standard xyz coordinates file
//In this version it cannot distinguish different elemental species
//cut_md is written in c++17.
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <omp.h>


using namespace std;

//CONSTANTS
//Set mass of the atoms  in Kg
float mass_Ar {39.948*1.66053*pow(10,-27)};
//Set K di boltzmann 1.380649×10−23 J/K = m2⋅kg/(s2⋅K)
float k_boltzmann {1.380649*pow(10, -23)};

/*                           parameters of the LJ potential                          */
float eta {1.0104}; // eta parameter, in eV
float sigma {2.0}; // sigma parameter, in Ang
float cut_off_distance {9.0};// cut_off distance for computing forces and energy


//classes
class Atom {
    public:

    //attributes
    string element_name;
    vector <float> position {0,0,0};
    vector <float> force {0,0,0};
    vector <float> velocity {0,0,0};


    // methods

    // LJ_energy: energy between atom pair
    // input atom object, output LJ energy as a float
    float LJ_energy (Atom atom_pair){
        float energy_lennard_jones {0};
        float square_dist {0};
        float d {0};
        for (unsigned i {0}; i<=2;++i)
            {
            square_dist+=pow(position.at(i)-atom_pair.position.at(i),2);
            }
        d = sqrt(square_dist);
        energy_lennard_jones = eta*(pow((sigma/d),12)-pow((sigma/d),6));
        return energy_lennard_jones;
    };

    // calculation of the absolute value of the distance, d
    float distance (Atom atom_pair){
        float square_dist {0};
        float tmp_xyz {0};
        float d {0};
        vector <float> xyz_dist {0,0,0};
        for (unsigned xyz {0}; xyz<=2;++xyz){
            tmp_xyz=atom_pair.position.at(xyz)-position.at(xyz);
            square_dist+=pow(tmp_xyz,2);
        }
        d = sqrt(square_dist);
        return d;
    }

    // atom_atom_force_LJ: LJ force between an atom pair, in its x y z coordinates
    // input: atom object, output: vector with x y z components of the force
    vector <float> atom_atom_force_LJ (Atom atom_pair){
        float force_lennard_jones_tot {1};
        // calculation of the absolute value of the distance, d
        float square_dist {0};
        float tmp_xyz {0};
        float d {0};
        vector <float> xyz_dist {0,0,0};
        for (unsigned xyz {0}; xyz<=2;++xyz)
            {
            tmp_xyz=atom_pair.position.at(xyz)-position.at(xyz);
            square_dist+=pow(tmp_xyz,2);
            xyz_dist.at(xyz)=tmp_xyz;
            }
        d = sqrt(square_dist);
        //absolute value of the LJ force
        force_lennard_jones_tot = eta * (-12*pow(sigma,12)/pow(d,13)+6*pow(sigma,6)/pow(d,7));
        //redefinition of the atom_pair coordinate considering
        // the atom under esamination as the center 0 0 0.
        // the projection on x y z axes of LJ force is proportional to those values of x y z
        vector <float> output {};
        float x;
        float tmp;
        for (int xyz {0}; xyz<3; xyz++)
            {
            tmp = xyz_dist.at(xyz)*force_lennard_jones_tot/d;
            output.push_back(tmp);
            }
        return output;
        };
};

//PROTOTIPES FUNCTIONS
vector <Atom> import_atoms (string path_input);
vector <vector <unsigned>> screen_couples_interaction (vector <Atom> atoms_vector_list);
vector <Atom> compute_force (vector <Atom> atoms_vector_list, vector <vector <unsigned>> couples);
vector <Atom> starting_velocities (vector <Atom> atoms_vector_list, float T);
vector <Atom> compute_positions (vector <Atom> atoms_vector_list, float time_step);
vector <Atom> compute_velocity_and_forces (vector <Atom> atoms_vector_list,
float time_step, vector <vector <unsigned>> couples);
vector <Atom> termostat (vector <Atom> atoms_vector_list, float T);
float compute_energy (vector <Atom> atoms_vector_list, vector <vector <unsigned>> couples);
float compute_temperature (vector <Atom> atoms_vector_list);
void print_info_toscreen (vector <Atom> atoms_vector_list, int i, vector <vector <unsigned>> couples);
void print_info_tofile (vector <Atom> atoms_vector_list, ofstream& output_file);




int main (int argc, char *argv[] ) {
//initial time
chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

/*                           parameters of the dynamics                          */
//units of measurements: space: Ang. time: fs, energy: eV
//Set Temperature in K
float temperature {300.0};
//time step in fs
float time_step{1};
//MD steps
int steps_temostat {20},  // how often switch on the ternostat
steps_couples_calc {20}, // how often to recalculate the atoms pairwise interaciotn
steps_dynamic {10000}; // how many steps of MD?
// input file name
string input_file_name {"test.txt"};
// output file name
string output_file_name {"output.xyz"};



//define the vector of Atom objects
vector <Atom> atoms_vector_list;
// read and extract the data of a input file
atoms_vector_list = import_atoms(input_file_name);
//initialize random velocities
atoms_vector_list = starting_velocities(atoms_vector_list, temperature);
//open output file
ofstream output_file;
output_file.open (output_file_name);
//compute the couples of atoms within the cutoff
vector <vector <unsigned>> couples;
couples=screen_couples_interaction(atoms_vector_list);
//compute initial forces
atoms_vector_list = compute_force(atoms_vector_list, couples);
    
/*                           printing options                          */
print_info_tofile (atoms_vector_list, output_file);
print_info_toscreen(atoms_vector_list,0, couples);

/////RUN THE DYNAMICs////////
for (int i {1};i <= steps_dynamic; i++)
{
    //compute positions
    atoms_vector_list = compute_positions(atoms_vector_list, time_step);
    //compute the new velocity and new forces
    atoms_vector_list = compute_velocity_and_forces(atoms_vector_list, time_step, couples);
    //every now and then termostat the MD!
    if (i % steps_temostat == 0)
        atoms_vector_list=termostat(atoms_vector_list, temperature);
    //every now and then recopute the couples!
    if (i % steps_couples_calc == 0)
        couples=screen_couples_interaction(atoms_vector_list);
    //printing options
    print_info_toscreen(atoms_vector_list,i,couples);
    print_info_tofile (atoms_vector_list, output_file);
chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
cout << "Wall time MD cycle= " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << endl;
}
//final time

return 0;
}


//FUNCTIONs

// import atoms from file
// INPUT: path to file
//OUTPUT : vector of Atom objects
vector <Atom> import_atoms (string path_input) {
    vector <Atom> atoms_vector_list;
    ifstream input_file;
    input_file.open(path_input);
    string line;
    while (getline (input_file,  line)){
        string line2;
        stringstream X(line);
        string a;
        float b, c, d;
        while (X >> a >> b >> c >> d)
            {
                Atom new_atom;
                new_atom.element_name=a;
                new_atom.position.at(0)=b;
                new_atom.position.at(1)=c;
                new_atom.position.at(2)=d;
                atoms_vector_list.push_back(new_atom);
            }
    }
    return atoms_vector_list;
}

//compute the distance and store the atoms pair which are within a certain distance cutoff
vector <vector <unsigned>> screen_couples_interaction (vector <Atom> atoms_vector_list){
    //create a list of atom-atom interactions to study within the cutoff
    vector <vector <unsigned>> couples;
    for (unsigned i {0}; i<atoms_vector_list.size(); i++)
        for (unsigned  j{i+1}; j<atoms_vector_list.size(); j++)
            if (atoms_vector_list.at(i).distance(atoms_vector_list.at(j)) < cut_off_distance)
                couples.push_back({i,j});
    return couples;
}
// compute and refresh the forces of each atom object
// INPUT: atom list, vector with couples of interacting atoms
// OUTPUT: atom list with refreshed forces values
vector <Atom> compute_force (vector <Atom> atoms_vector_list, vector <vector <unsigned>> couples) {
    //force matrix -- no diagonal terms only one side due to symmetry
    //"array" full of zeros
    vector <float> force (3);
    vector <vector <float>> force_vector (couples.size(),force);
    //loop on the list of couples.
    #pragma omp parallel  for
            for (unsigned k = 0; k<couples.size(); k++){
                force_vector.at(k)=
                atoms_vector_list.at(couples.at(k).at(0)).
                atom_atom_force_LJ(atoms_vector_list.at(couples.at(k).at(1)));
            }
    //clean forces
    vector <float> v_null {0,0,0};
    for (int i {0}; i < atoms_vector_list.size(); i++)
        atoms_vector_list.at(i).force=v_null;
    // sum over colum and row
    #pragma omp parallel shared (atoms_vector_list)
    {
        #pragma omp for
        for (unsigned k=0; k<couples.size(); k++)
                for (unsigned xyz {0}; xyz <3; xyz++){
                        atoms_vector_list.at(couples.at(k).at(0)).force.at(xyz)+=
                    force_vector.at(k).at(xyz);
                        atoms_vector_list.at(couples.at(k).at(1)).force.at(xyz)-=
                    force_vector.at(k).at(xyz);}
    }
return atoms_vector_list;
}

// set random initial velocities to the atoms coherent to desider T
// INPUT: atom list with v=0
//OUTPUT: atom list with refreshed velocities
vector <Atom> starting_velocities (vector <Atom> atoms_vector_list, float T){
    float max_velocity_in_single_direction {0};
    //max velocity considering the T = 300 K in m/s and converted to  A/fs  (*scaling factor 10^-5)
    max_velocity_in_single_direction = pow(3*T*k_boltzmann/mass_Ar,0.5)*pow(10,-5);
    //generate a truple of number in which one is +1 or -1 and the other two are zeros
    vector <vector <float>> ntuple {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
    int random;
// srand((unsigned) time(NULL));
    for (unsigned i {0} ; i<atoms_vector_list.size(); i++){
        random = rand() %  5;
        for (unsigned  xyz{0} ; xyz<=2; xyz++){
            atoms_vector_list.at(i).velocity.at(xyz) = ntuple.at(random).at(xyz)*max_velocity_in_single_direction;
        }
    }
    return atoms_vector_list;
}

//MOVE: function which uses velet intergration to move the atoms of one step
//INPUT: vector of atoms, time step
//OUTPUT:vector of atoms with positions uploaded
vector <Atom> compute_positions (vector <Atom> atoms_vector_list, float time_step){
    for (unsigned i {0} ; i<atoms_vector_list.size(); i++){
        for (unsigned  xyz{0} ; xyz<=2; xyz++){
            //convert the force from eV/A to A*Kg/fs^2
                float force_coverted {atoms_vector_list.at(i).force.at(xyz)*1.602176634*pow(10,-29)};
            atoms_vector_list.at(i).position.at(xyz)+=atoms_vector_list.at(i).velocity.at(xyz)*time_step+0.5*pow(time_step,2)*force_coverted/mass_Ar;
        }
    }
    return atoms_vector_list;
    }

//VELOCITY VERLET FOR FORCES AND VELOCITIES
vector <Atom> compute_velocity_and_forces (vector <Atom> atoms_vector_list, float time_step,  vector <vector <unsigned>> couples){
    //saving current force on a vector
    vector <vector <float>> current_forces;
    for (unsigned i {0} ; i<atoms_vector_list.size(); i++)
        current_forces.push_back(atoms_vector_list.at(i).force);
    //recomputing forces
    atoms_vector_list = compute_force(atoms_vector_list,couples);
    //computing velocity using old and new forces
    for (unsigned i {0} ; i<atoms_vector_list.size(); i++){
        for (unsigned  xyz{0} ; xyz<=2; xyz++){
            //convert the force from eV/A to A*Kg/fs^2
            float force_coverted_new {atoms_vector_list.at(i).force.at(xyz)*1.602176634*pow(10,-29)};
            float force_coverted_old {current_forces.at(i).at(xyz)*1.602176634*pow(10,-29)};
            atoms_vector_list.at(i).velocity.at(xyz)=atoms_vector_list.at(i).velocity.at(xyz)+(force_coverted_new+force_coverted_old)/mass_Ar*time_step*0.5;
        }
    }
    return atoms_vector_list;
}

//COMPUTE  the temperature of the simulation
float compute_temperature (vector <Atom> atoms_vector_list){
    float sum_atoms_velocity {0};
    for (unsigned i {0} ; i<atoms_vector_list.size(); i++){
        float tmp {0};
        for (unsigned  xyz{0} ; xyz<=2; xyz++){
            tmp+=pow(atoms_vector_list.at(i).velocity.at(xyz),2);
        }
        sum_atoms_velocity += tmp*mass_Ar*pow(10,10)/(3*k_boltzmann);
    }
    //return the T value converting the velocity with the correct unit of mesurement
    //and averaging the overall value for the number of atoms
    return sum_atoms_velocity/atoms_vector_list.size();
}

//RESCALE temperature
vector <Atom> termostat (vector <Atom> atoms_vector_list, float T){
     float scaling_factor {T/compute_temperature(atoms_vector_list)};
     for (unsigned i {0} ; i<atoms_vector_list.size(); i++){
        for (unsigned  xyz{0} ; xyz<=2; xyz++){
            atoms_vector_list.at(i).velocity.at(xyz)*= pow(scaling_factor,0.5);
        }
     }
     return atoms_vector_list;
}

 //COMPUTE total energy
float compute_energy (vector <Atom> atoms_vector_list, vector <vector <unsigned>> couples){
    float total_energy{0};
    #pragma omp parallel for reduction(+:total_energy)
                for (unsigned k=0; k<couples.size(); k++)
                        total_energy+=atoms_vector_list.at(couples.at(k).at(0)).LJ_energy(atoms_vector_list.at(couples.at(k).at(1)));
    return total_energy;
 }

//PRINTING SECTION
void print_info_toscreen (vector <Atom> atoms_vector_list, int i, vector <vector <unsigned>> couples){
    cout << "step " << i ;
    cout << " temperature "<< compute_temperature(atoms_vector_list);
    cout << " E= "<<compute_energy(atoms_vector_list,couples)<<" eV";
//    cout << endl;
//    for (unsigned i {0} ; i<atoms_vector_list.size(); i++){
//        cout << atoms_vector_list.at(i).element_name<<" ";
//        for (unsigned  xyz{0} ; xyz<=2; xyz++)
//            cout << atoms_vector_list.at(i).position.at(xyz) << " ";
//        cout <<" velocities ";
//        for (unsigned  xyz{0} ; xyz<=2; xyz++)
//            cout << atoms_vector_list.at(i).velocity.at(xyz) << " ";
//        cout <<" forces ";
//        for (unsigned  xyz{0} ; xyz<=2; xyz++)
//            cout << atoms_vector_list.at(i).force.at(xyz) << " ";
//        cout << endl;
//    }
//    cout << endl;
}

void print_info_tofile (vector <Atom> atoms_vector_list, ofstream& output_file){
output_file << atoms_vector_list.size()<<endl;
output_file <<"comment"<<endl;
    for (unsigned i {0} ; i<atoms_vector_list.size(); i++){
        output_file << atoms_vector_list.at(i).element_name<<" ";
        for (unsigned  xyz{0} ; xyz<=2; xyz++)
            output_file << atoms_vector_list.at(i).position.at(xyz) << " ";
//        output_file <<" velocities ";
//        for (unsigned  xyz{0} ; xyz<=2; xyz++)
//            output_file << atoms_vector_list.at(i).velocity.at(xyz) << " ";
//        output_file <<" forces ";
//        for (unsigned  xyz{0} ; xyz<=2; xyz++)
//            output_file << atoms_vector_list.at(i).force.at(xyz) << " ";
        output_file << endl;
    }
    cout << endl;
}

