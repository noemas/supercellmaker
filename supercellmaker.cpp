//__________________________________INPUT/OUTPUT__________________________________//
//Name of the input files/cells in quotation marks.
//Need to be in .xyz format or .vasp (POSCAR) format.
//If you want to rotate/shift a single file set the second input as an empty string ("").
#define INPUT1 ""
#define INPUT2 ""

//Name of the output file and output format ("POSCAR" or "xyz") in quotation marks.
#define OUTPUT ""
#define OUTPUT_FORMAT "POSCAR"
//_________________________________________________________________________________//


//_________________________________CELL ROTATIONS__________________________________//
//Define the rotation angles for the two cells.
//These are not Euler angles but the rotation angles around global rotation vectors defined further below.
//Rotations are made before translations.
//Disabled if all angles are 0.

//Rotation of 1st cell (in degrees)
#define ALPHA1 0
#define BETA1  0
#define GAMMA1 0
//Rotation of 2nd cell (in degrees)
#define ALPHA2 0
#define BETA2 0
#define GAMMA2 0

//Define custom rotation axes (just the direction is relevant, i.e. length doesn't matter) in cartesian coordinates.
//Rotation angles ALPHA1/2 act on on ROT_A, BETA1/2 on ROT_B, GAMMA1/2 on ROT_C 
#define ROT_A_X 1
#define ROT_A_Y 0
#define ROT_A_Z 0

#define ROT_B_X 0
#define ROT_B_Y 1
#define ROT_B_Z 0

#define ROT_C_X 0
#define ROT_C_Y 0
#define ROT_C_Z 1
//_________________________________________________________________________________//

//________________________________CELL TRANSLATIONS________________________________//
//Cell translations in Angstrom and cartesian coordinates. Translations are made after the rotations.
//Disabled if all translations are 0.

//Translation of 1st cell
#define SHIFT_X1 0
#define SHIFT_Y1 0
#define SHIFT_Z1 0

//Translation of 2nd cell
#define SHIFT_X2 0
#define SHIFT_Y2 0
#define SHIFT_Z2 0
//_________________________________________________________________________________//



//______________________________UNIT CELL TO BE CUT________________________________//
//Define a unit cell here that will be cut out of the combined input files. Cutting will only take place if WALL_UC = false
//(i.e. if you don't generate a domain wall). In addition to the final file, a file called cut.xyz will be generated after the cut
//where the corners of the unit cells are denoted by additional atoms of element (Uc -> unit cell) and all deleted atoms are set to 
//a new element (De -> deleted).

//Unit cell lattice vectors in Angstrom and cartesian coordinates (just the direction is relevant, i.e. length doesn't matter)
#define A_X 1
#define A_Y 0
#define A_Z 0

#define B_X 0
#define B_Y 1
#define B_Z 0

#define C_X 0
#define C_Y 0
#define C_Z 1

//Lattice parameters in Angstrom.
//Cutting in a direction is disabled if the lattice parameter in this direction is 0.
#define A 0
#define B 0
#define C 0

//Translations in cartesian coordinates necessary to shift the origin of the unit cell to the true origin (0, 0, 0).
#define O_X 0
#define O_Y 0
#define O_Z 0

//Automatic determination of the unit cell.
//This is very crude and only possible if:
//	1. the structure is at least orthorhombic
//	2. the atomic structure after merging contains only the atoms in the unit cell (i.e. no cuts necessary)
//	3. there is at least one atom at the 0 coordinate in every direction
#define AUTO_UC true
//_________________________________________________________________________________//

//_____________________________DOMAIN_WALL_GENERATION______________________________//
//This will only yield reasonable results if the combination of the two unit cells is reasonable, i.e. doesn't create a cell with more than 8 faces.
//set WALL_UC to true if you want to create a domain wall using the two input structures. If WALL_UC is set to true, unit cell cutting feature is disabled.
#define WALL_UC false
//direction of the domain wall given as the lattice parameter index, e.g. if WALL_DIR is 2 the domain wall will be constructed lattice vector b.
#define WALL_DIR 2
//sizes of the two domains. Syntax: "n1 n2 n3", where n1 is the number of unit cells along the first lattice parameter etc.
#define DOMAIN_SIZE1 "1 4 1"
#define DOMAIN_SIZE2 "1 4 1"

//____________________________________DOUBLES______________________________________//
//Double atoms tolerance in Angstrom. Any atoms of the same type with a distance lower than TOL will be considered as doubles (the same) and therefore be averaged.
//If atoms have a distance greater than TOL but can be joined by atoms with a distance lower than TOL, then all these atoms are considered as doubles.
//Also true across the cell boundaries.
//If you want to average two distinct atoms try to set TOL_MIN slightly below and TOL_MAX slightly above the bond length. 
//Disabled if TOL_MIN is larger than TOL_MAX or if they have the same value.
#define TOL_MIN 0
#define TOL_MAX 0.8
//_________________________________________________________________________________//


//____________________________________SORTING______________________________________//
//Sort the atoms along a direction by defining it in cartesian coordinates.
//Disabled if sorting direction is [0,0,0].
#define SORT_X 1
#define SORT_Y 0
#define SORT_Z 0

//Choose if the sorted atoms need to be in blocked form (atoms of the same types are sorted along specified direction but grouped in blocks of the same type).
//BLOCK needs to be true if output format is POSCAR.
#define BLOCK true

//______________________________________END________________________________________//
//______________________________NOW_COMPILE_AND_RUN________________________________//
//_________________________________________________________________________________//

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<algorithm>
#include<cmath>
#include<sstream>
#include<functional>
#include<thread>
#include<mutex>


//_________________________________________________________________________________________________________//

const std::string input_file_1 = (std::string)INPUT1;
const std::string input_file_2 = (std::string)INPUT2;
const std::string output_file = (std::string)OUTPUT;
const std::string output_format = (std::string)OUTPUT_FORMAT;
const double shift_x_1 = (double)SHIFT_X1;
const double shift_y_1 = (double)SHIFT_Y1;
const double shift_z_1 = (double)SHIFT_Z1;
const double shift_x_2 = (double)SHIFT_X2;
const double shift_y_2 = (double)SHIFT_Y2;
const double shift_z_2 = (double)SHIFT_Z2;
const double alpha1 = (double)ALPHA1;
const double beta1 = (double)BETA1;
const double gamma1 = (double)GAMMA1;
const double alpha2 = (double)ALPHA2;
const double beta2 = (double)BETA2;
const double gamma2 = (double)GAMMA2;
const double rot_a_x = (double)ROT_A_X;
const double rot_a_y = (double)ROT_A_Y;
const double rot_a_z = (double)ROT_A_Z;
const double rot_b_x = (double)ROT_B_X;
const double rot_b_y = (double)ROT_B_Y;
const double rot_b_z = (double)ROT_B_Z;
const double rot_c_x = (double)ROT_C_X;
const double rot_c_y = (double)ROT_C_Y;
const double rot_c_z = (double)ROT_C_Z;
const double tol_min = (double)TOL_MIN;
const double tol_max = (double)TOL_MAX;
const double sort_x = (char)SORT_X;
const double sort_y = (char)SORT_Y;
const double sort_z = (char)SORT_Z;
const bool block = (bool)BLOCK;
const double a_x = (double)A_X;
const double a_y = (double)A_Y;
const double a_z = (double)A_Z;
const double b_x = (double)B_X;
const double b_y = (double)B_Y;
const double b_z = (double)B_Z;
const double c_x = (double)C_X;
const double c_y = (double)C_Y;
const double c_z = (double)C_Z;
const double a = (double)A;
const double b = (double)B;
const double c = (double)C;
const double o_x = (double)O_X;
const double o_y = (double)O_Y;
const double o_z = (double)O_Z;
const bool auto_uc = (bool)AUTO_UC;
const bool wall_uc = (bool)WALL_UC;
const int wall_dir = (int)WALL_DIR;
const std::string domain_size_1 = (std::string) DOMAIN_SIZE1;
const std::string domain_size_2 = (std::string) DOMAIN_SIZE2;
const bool open_window = true;
void keep_window_open(bool _open_window)
{
	if (_open_window)
	{
		std::cout << "Press enter to exit" << std::endl;
		std::cin.get();
	}
}

//_________________________________________________________________________________________________________//


//simple 3d vector struct for basic vector arithmetic
struct vector3D
{
	vector3D();
	vector3D(double, double, double);
	vector3D(double, double, double, double);
	double x, y, z, length;
	void calc_length();
	double calc_dist(const vector3D&) const;
	double calc_angle(const vector3D&) const;
	vector3D operator*(double) const;
	double operator*(const vector3D&) const;
	vector3D operator+(const vector3D&) const;
	vector3D operator-(const vector3D&) const;
};

class atomic_structure
{
public:
	atomic_structure();
	atomic_structure(const atomic_structure&);
	atomic_structure(const std::string, const std::string);
	atomic_structure& operator=(const atomic_structure&);
	void read_xyz();
	void read_poscar();
	void rotate(double, double, double, const vector3D&, const vector3D&, const vector3D&);
	void shift(double, double, double);
	atomic_structure merge(const atomic_structure&);
	void cut(double, double, double);
	void write_xyz(std::string);
	void write_poscar(std::string);
	void average_atoms(double, double, std::string);
	void sort(const vector3D&, bool);
	void set_lat(const vector3D&, const vector3D&, const vector3D&);
	void find_uc();
	std::vector<double> min_max();
	void generate_wall(const atomic_structure&, std::string, std::string, int, double, double);
	void complete_uc(const std::vector<double>&, int, double, double);
	std::vector<double> calc_T_inv();
private:
	std::ifstream ifs;	//input file stream
	int num_of_atoms;	//total number of atoms in the structure
	std::vector<std::string> elements;	//vector containing the present element types							
	std::vector<std::string> elements_vec;	//vector containing the element types of the individual atoms
	std::vector<std::vector<double>> position_matrix;	//position matrix containing the cartesian atom coordinates in angstrom			
	vector3D lat_a, lat_b, lat_c;	//lattice vectors in angstrom				   
	std::vector<int> find_neighbours(int, double, double);
	std::vector<int> find_neighbours_mic(int, double, double);
	std::vector<int> get_neighbour_web(std::vector<int>, double, double, std::string);
	void move_images(int, const std::vector<int>&);
};

//_________________________________________________________________________________________________________//

vector3D::vector3D()
	:x(0), y(0), z(0), length(0)
{}

vector3D::vector3D(double _x, double _y, double _z)
	: x(_x), y(_y), z(_z)
{
	this->calc_length();
}

//construct vector of specified length
vector3D::vector3D(double _x, double _y, double _z, double _length)
	: x(_x), y(_y), z(_z)
{
	this->calc_length();
	double scale = _length / this->length;
	this->x *= scale;
	this->y *= scale;
	this->z *= scale;
	this->calc_length();
}

//calculates and returns the length/norm of the vector
void vector3D::calc_length()
{
	this->length = sqrt(pow(this->x, 2) + pow(this->y, 2) + pow(this->z, 2));
}

//calculates and returns the distance between two vectors/points
double vector3D::calc_dist(const vector3D& _vec2) const
{
	return sqrt(pow(this->x - _vec2.x, 2) + pow(this->y - _vec2.y, 2) + pow(this->z - _vec2.z, 2));
}

//calculates and returns the angle between two vectors in radians
double vector3D::calc_angle(const vector3D& _vec2) const
{
	return acos(*this * _vec2 / (this->length * _vec2.length));
}

//basic arithmetic operators
vector3D vector3D::operator*(double _scalar) const
{
	vector3D new_vec = *this;
	new_vec.x *= _scalar;
	new_vec.y *= _scalar;
	new_vec.z *= _scalar;
	new_vec.calc_length();
	return new_vec;
}

double vector3D:: operator*(const vector3D& _vec2) const
{
	return this->x*_vec2.x + this->y*_vec2.y + this->z*_vec2.z;
}

vector3D vector3D::operator+(const vector3D& _vec2) const
{
	vector3D new_vec = *this;
	new_vec.x += _vec2.x;
	new_vec.y += _vec2.y;
	new_vec.z += _vec2.z;
	new_vec.calc_length();
	return new_vec;
}

vector3D vector3D::operator-(const vector3D& _vec2) const
{
	vector3D new_vec = *this;
	new_vec.x -= _vec2.x;
	new_vec.y -= _vec2.y;
	new_vec.z -= _vec2.z;
	new_vec.calc_length();
	return new_vec;
}


atomic_structure::atomic_structure()
	:num_of_atoms(0), elements(std::vector<std::string>(0)), elements_vec(std::vector<std::string>(0)),
	position_matrix(std::vector<std::vector<double>>(0, std::vector<double>(3, 0))),
	lat_a(0, 0, 0), lat_b(0, 0, 0), lat_c(0, 0, 0)
{}

atomic_structure::atomic_structure(const atomic_structure& second_file)
	: num_of_atoms(second_file.num_of_atoms), elements(second_file.elements), elements_vec(second_file.elements_vec),
	position_matrix(second_file.position_matrix), lat_a(second_file.lat_a), lat_b(second_file.lat_b), lat_c(second_file.lat_c)
{}

//constructor for atomic structure which opens and reads a file. If file extension is not recognized, an attempt to read
//the file in POSCAR format will be made
atomic_structure::atomic_structure(std::string _file, std::string _input_format)
	:ifs(_file)
{
	if (ifs.good())
	{
		if (_input_format == "xyz")
		{
			this->read_xyz();
		}
		else
		{
			this->read_poscar();
		}
	}
	else if (!ifs.good() && _file.empty())
	{
		std::string error_msg("Error: Empty file.");
		throw(error_msg);
	}
}

atomic_structure& atomic_structure::operator=(const atomic_structure& _second_file)
{
	this->num_of_atoms = _second_file.num_of_atoms;
	this->elements = _second_file.elements;
	this->elements_vec = _second_file.elements_vec;
	this->position_matrix = _second_file.position_matrix;
	this->lat_a = _second_file.lat_a;
	this->lat_b = _second_file.lat_b;
	this->lat_c = _second_file.lat_c;
	return *this;
}

void atomic_structure::read_xyz()
{
	this->ifs >> this->num_of_atoms; //get number of atoms
	std::string eol;
	getline(this->ifs, eol);		//get end of 1st line
	getline(this->ifs, eol);		//discard comment (2nd) line

	this->elements_vec.resize(this->num_of_atoms);
	this->position_matrix.resize(this->num_of_atoms, std::vector<double>(3, 0));
	std::string el;
	double x_pos(0), y_pos(0), z_pos(0);
	int index(0);
	//read atom type and position data
	while (this->ifs >> el >> x_pos >> y_pos >> z_pos)
	{
		this->elements_vec[index] = el;
		this->position_matrix[index][0] = x_pos;
		this->position_matrix[index][1] = y_pos;
		this->position_matrix[index][2] = z_pos;
		if (std::find(this->elements.begin(), this->elements.end(), el) == this->elements.end())
		{
			this->elements.push_back(el);
		}
		index++;
		if (index > this->num_of_atoms)
		{
			std::string error_msg("Error: Number of atoms in input file does not match number of atom positions.");
			throw(error_msg);
		}
	}

	this->ifs.close();
	if (index != this->num_of_atoms)
	{
		std::string error_msg("Error: Number of atoms in input file does not match number of atom positions.");
		throw(error_msg);
	}
	if (this->num_of_atoms == 0 || index == 0)
	{
		std::string error_msg("Error: Input file either contains no atoms or format is not recognized.");
		throw(error_msg);
	}
}

void atomic_structure::read_poscar()
{
	//discard comment in 1st line
	std::string comment;
	getline(this->ifs, comment);
	//get lattice scaling factor
	double scaling(0);
	this->ifs >> scaling;
	//get lattice vectors
	double a_x(0), a_y(0), a_z(0), b_x(0), b_y(0), b_z(0), c_x(0), c_y(0), c_z(0);
	this->ifs >> a_x >> a_y >> a_z >> b_x >> b_y >> b_z >> c_x >> c_y >> c_z;
	vector3D a(a_x, a_y, a_z);
	vector3D b(b_x, b_y, b_z);
	vector3D c(c_x, c_y, c_z);
	this->lat_a = a*scaling;
	this->lat_b = b*scaling;
	this->lat_c = c*scaling;

	//get next line (either atom types or number of atoms)
	std::string line;
	getline(this->ifs, line); //end of line
	getline(this->ifs, line);
	std::stringstream ss(line);
	int num_of_atoms_el(0);
	ss >> num_of_atoms_el;
	//if this fails the line contains element types
	if (ss.fail())
	{
		//get the elements
		ss.clear();
		std::string el;
		while (ss >> el)
		{
			this->elements.push_back(el);
		}
		//get the next line containing the number of atoms
		getline(this->ifs, line);
		ss.str(line);
	}
	else { ss.str(line); }
	ss.clear();

	if (this->elements.size() == 0)
	{
		std::cout << "Warning: No element types specified in POSCAR" << std::endl;
	}
	//read the number of atoms per element and set the elements array accordingly
	int index(0);
	while (ss >> num_of_atoms_el)
	{
		this->num_of_atoms += num_of_atoms_el;
		//give the elements placeholder names if no elements were read
		if (this->elements.size() == index)
		{
			this->elements.push_back("E" + std::to_string(index));
		}
		for (int i = 0; i < num_of_atoms_el; i++)
		{
			this->elements_vec.push_back(this->elements[index]);
		}
		index++;
	}

	this->position_matrix.resize(this->num_of_atoms, std::vector<double>(3, 0));

	//get the type of coordinate
	bool direct(true);
	this->ifs >> line;

	switch (line[0])
	{
	case 'S':
	{
		std::string error_msg("Error: Selective Dynamics not supported.");
		throw(error_msg);
		break;
	}
	case 'C':
	{
		direct = false;
	}
	break;
	case 'D':
		break;
	default:
	{
		std::string error_msg("Error: Expecting Cartesian or Direct coordinates specification on line 8.");
		throw(error_msg);
		break;
	}
	}

	//read the positions
	double x_pos(0), y_pos(0), z_pos(0);
	index = 0;
	//read atom type and position data
	for (int i = 0; i < this->num_of_atoms; i++)
	{
		this->ifs >> x_pos >> y_pos >> z_pos;
		vector3D atom_pos;
		if (direct)
		{
			atom_pos = this->lat_a * x_pos + this->lat_b * y_pos + this->lat_c * z_pos;
		}
		else
		{
			atom_pos.x = x_pos;
			atom_pos.y = y_pos;
			atom_pos.z = z_pos;
		}
		this->position_matrix[i] = { atom_pos.x, atom_pos.y, atom_pos.z };
		index++;
	}
}

//rotate the structure around the specified roation axes by the specified amount of degrees in degrees
void atomic_structure::rotate(double _alpha, double _beta, double _gamma,
	const vector3D& _rot1, const vector3D& _rot2, const vector3D& _rot3)
{
	const double pi = 4 * atan(1);
	const double rot_x = 2 * pi / 360 * _alpha;
	const double rot_y = 2 * pi / 360 * _beta;
	const double rot_z = 2 * pi / 360 * _gamma;
	std::vector<double> angles{ rot_x, rot_y, rot_z };
	std::vector<vector3D> lat_vecs(3);
	lat_vecs[0] = _rot1;
	lat_vecs[1] = _rot2;
	lat_vecs[2] = _rot3;
	std::vector<double> ROT(9, 0);

	for (int i = 2; i >= 0; i--)
	{
		double angle = angles[i];
		vector3D lat_vec = lat_vecs[i];
		ROT[0] = cos(angle) + pow(lat_vec.x, 2)*(1 - cos(angle));
		ROT[1] = lat_vec.x * lat_vec.y*(1 - cos(angle)) - lat_vec.z*sin(angle);
		ROT[2] = lat_vec.x * lat_vec.z*(1 - cos(angle)) + lat_vec.y *sin(angle);
		ROT[3] = lat_vec.x * lat_vec.y * (1 - cos(angle)) + lat_vec.z * sin(angle);
		ROT[4] = cos(angle) + pow(lat_vec.y, 2)*(1 - cos(angle));
		ROT[5] = lat_vec.y * lat_vec.z*(1 - cos(angle)) - lat_vec.x * sin(angle);
		ROT[6] = lat_vec.x*lat_vec.z*(1 - cos(angle)) - lat_vec.y*sin(angle);
		ROT[7] = lat_vec.y* lat_vec.z*(1 - cos(angle)) + lat_vec.x*sin(angle);
		ROT[8] = cos(angle) + pow(lat_vec.z, 2)*(1 - cos(angle));
		//rotate the atomic positions
		for (int j = 0; j < this->num_of_atoms; j++)
		{
			double old_x = this->position_matrix[j][0];
			double old_y = this->position_matrix[j][1];
			double old_z = this->position_matrix[j][2];

			this->position_matrix[j][0] = ROT[0] * old_x + ROT[1] * old_y + ROT[2] * old_z;
			this->position_matrix[j][1] = ROT[3] * old_x + ROT[4] * old_y + ROT[5] * old_z;
			this->position_matrix[j][2] = ROT[6] * old_x + ROT[7] * old_y + ROT[8] * old_z;
		}
		//rotate the lattice vectors
		double old_a_x = this->lat_a.x;
		double old_a_y = this->lat_a.y;
		double old_a_z = this->lat_a.z;
		double old_b_x = this->lat_b.x;
		double old_b_y = this->lat_b.y;
		double old_b_z = this->lat_b.z;
		double old_c_x = this->lat_c.x;
		double old_c_y = this->lat_c.y;
		double old_c_z = this->lat_c.z;
		this->lat_a.x = ROT[0] * old_a_x + ROT[1] * old_a_y + ROT[2] * old_a_z;
		this->lat_a.y = ROT[3] * old_a_x + ROT[4] * old_a_y + ROT[5] * old_a_z;
		this->lat_a.z = ROT[6] * old_a_x + ROT[7] * old_a_y + ROT[8] * old_a_z;
		this->lat_b.x = ROT[0] * old_b_x + ROT[1] * old_b_y + ROT[2] * old_b_z;
		this->lat_b.y = ROT[3] * old_b_x + ROT[4] * old_b_y + ROT[5] * old_b_z;
		this->lat_b.z = ROT[6] * old_b_x + ROT[7] * old_b_y + ROT[8] * old_b_z;
		this->lat_c.x = ROT[0] * old_c_x + ROT[1] * old_c_y + ROT[2] * old_c_z;
		this->lat_c.y = ROT[3] * old_c_x + ROT[4] * old_c_y + ROT[5] * old_c_z;
		this->lat_c.z = ROT[6] * old_c_x + ROT[7] * old_c_y + ROT[8] * old_c_z;
	}
}

//shift the structure by the specified translation in angstrom
void atomic_structure::shift(double _a, double _b, double _c)
{
	for (int i = 0; i < this->num_of_atoms; i++)
	{
		this->position_matrix[i][0] += _a;
		this->position_matrix[i][1] += _b;
		this->position_matrix[i][2] += _c;
	}
}

//merge two atomic structures. The lattice parameters are lost if present.
atomic_structure atomic_structure::merge(const atomic_structure& _second_file)
{
	atomic_structure merged_file;
	merged_file.num_of_atoms = this->num_of_atoms + _second_file.num_of_atoms;
	merged_file.elements_vec = this->elements_vec;
	merged_file.position_matrix = this->position_matrix;
	for (int i = 0; i < _second_file.num_of_atoms; i++)
	{
		std::string el = _second_file.elements_vec[i];
		std::vector<double> pos = _second_file.position_matrix[i];
		int index = std::find(merged_file.elements_vec.begin(), merged_file.elements_vec.end(), el) - merged_file.elements_vec.begin();
		if (index == merged_file.elements_vec.size())
		{
			merged_file.elements_vec.insert(merged_file.elements_vec.end(), el);
			merged_file.position_matrix.insert(merged_file.position_matrix.end(), pos);
		}
		else
		{
			merged_file.elements_vec.insert(merged_file.elements_vec.begin() + index, el);
			merged_file.position_matrix.insert(merged_file.position_matrix.begin() + index, pos);
		}
	}
	merged_file.elements = this->elements;
	for (unsigned int i = 0; i < _second_file.elements.size(); i++)
	{
		if (std::find(merged_file.elements.begin(), merged_file.elements.end(), _second_file.elements[i]) == merged_file.elements.end())
		{
			merged_file.elements.insert(merged_file.elements.end(), _second_file.elements[i]);
		}
	}
	return merged_file;
}

//cut out out all atoms that lie in the specified lattice parameters of the atomic structure. The input is the necessary translation which
//shifts the atoms that are to be cut in the unit cell specified by the lattice parameters
void atomic_structure::cut(double _x_shift, double _y_shift, double _z_shift)
{
	atomic_structure cut_display = *this;

	//shift the atoms to the origin
	this->shift(_x_shift, _y_shift, _z_shift);

	//get the direct coordinates
	std::vector<double> T_inv = this->calc_T_inv();
	double x_direct(0), y_direct(0), z_direct(0);

	int i(0);
	int cut_atoms(0);
	while (i != this->num_of_atoms)
	{
		bool cut(false);
		x_direct = T_inv[0] * this->position_matrix[i][0] + T_inv[1] * this->position_matrix[i][1] + T_inv[2] * this->position_matrix[i][2];
		y_direct = T_inv[3] * this->position_matrix[i][0] + T_inv[4] * this->position_matrix[i][1] + T_inv[5] * this->position_matrix[i][2];
		z_direct = T_inv[6] * this->position_matrix[i][0] + T_inv[7] * this->position_matrix[i][1] + T_inv[8] * this->position_matrix[i][2];
		if (x_direct < 0 || x_direct > 1 || y_direct < 0 || y_direct > 1 || z_direct < 0 || z_direct > 1)
		{
			this->position_matrix.erase(this->position_matrix.begin() + i);
			this->elements_vec.erase(this->elements_vec.begin() + i);
			this->num_of_atoms--;
			cut = true;
			cut_atoms++;
			if (cut_atoms == 1)
			{
				cut_display.elements.push_back("DE");
			}
			cut_display.elements_vec[i + cut_atoms - 1] = "DE";
		}
		if (!cut)
		{
			i++;
		}
	}

	cut_display.elements.push_back("UC");
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			for (int k = 0; k < 2; k++)
			{
				cut_display.elements_vec.push_back("UC");
				vector3D corner = this->lat_a * i + this->lat_b * j + this->lat_c * k;
				corner.x -= _x_shift;
				corner.y -= _y_shift;
				corner.z -= _z_shift;
				cut_display.position_matrix.push_back({ corner.x, corner.y, corner.z });
				cut_display.num_of_atoms++;
			}
		}
	}

	cut_display.write_xyz("cut.xyz");

}

//write the atomic structure to an xyz file
void atomic_structure::write_xyz(std::string _file_name)
{
	std::ofstream ofs(_file_name);
	ofs << this->num_of_atoms << std::endl;
	for (unsigned int i = 0; i < this->elements.size(); i++)
	{
		ofs << this->elements[i] << '\t';
	}
	ofs << std::endl;
	for (int i = 0; i < this->num_of_atoms; i++)
	{
		ofs << this->elements_vec[i] << '\t' << this->position_matrix[i][0]
			<< '\t' << this->position_matrix[i][1] << '\t' << this->position_matrix[i][2] << std::endl;
	}

	ofs.close();
}

//write the atoms to a poscar file
void atomic_structure::write_poscar(std::string _file_name)
{
	std::vector<int> elements_num(this->elements.size(), 0);
	for (int i = 0; i < this->num_of_atoms; i++)
	{
		std::string el = this->elements_vec[i];
		int index = std::find(this->elements.begin(), this->elements.end(), el) - this->elements.begin();
		elements_num[index]++;
	}
	std::ofstream ofs(_file_name);

	for (unsigned int i = 0; i < this->elements.size(); i++)
	{
		ofs << elements[i] << " ";
	}

	ofs << std::endl << 1.0 << std::endl;

	ofs << this->lat_a.x << '\t' << this->lat_a.y << '\t' << this->lat_a.z << std::endl
		<< this->lat_b.x << '\t' << this->lat_b.y << '\t' << this->lat_b.z << std::endl
		<< this->lat_c.x << '\t' << this->lat_c.y << '\t' << this->lat_c.z << std::endl;

	for (unsigned int i = 0; i < this->elements.size(); i++)
	{
		ofs << elements[i] << " ";
	}

	ofs << std::endl;

	for (unsigned int i = 0; i < this->elements.size(); i++)
	{
		ofs << elements_num[i] << " ";
	}

	ofs << std::endl;
	ofs << "Direct" << std::endl;

	//get the direct coordinates and write them
	std::vector<double> T_inv = this->calc_T_inv();
	double x_direct(0), y_direct(0), z_direct(0);
	for (int i = 0; i < this->num_of_atoms; i++)
	{
		x_direct = T_inv[0] * this->position_matrix[i][0] + T_inv[1] * this->position_matrix[i][1] + T_inv[2] * this->position_matrix[i][2];
		y_direct = T_inv[3] * this->position_matrix[i][0] + T_inv[4] * this->position_matrix[i][1] + T_inv[5] * this->position_matrix[i][2];
		z_direct = T_inv[6] * this->position_matrix[i][0] + T_inv[7] * this->position_matrix[i][1] + T_inv[8] * this->position_matrix[i][2];
		ofs << '\t' << x_direct << '\t'
			<< y_direct << '\t'
			<< z_direct << std::endl;
	}
	ofs.close();
}

//average out atoms with distances which fall in the specified range. Last input is the output format.
void atomic_structure::average_atoms(double _tol_min, double _tol_max, std::string _format)
{
	int i = 0;
	bool bad_avg(false);
	while (i < this->num_of_atoms)
	{
		std::vector<int> neighbours;
		int num_of_nb(0);
		//find the neighbours of atom i
		if (_format == "xyz")
		{
			neighbours = this->find_neighbours(i, _tol_min, _tol_max);
		}
		else if (_format == "POSCAR")
		{
			neighbours = this->find_neighbours_mic(i, _tol_min, _tol_max);
		}
		num_of_nb = neighbours.size();
		//if there are no neighbours skip to the next atom
		if (num_of_nb == 0)
		{
			i++;
			continue;
		}
		//else find the spanning web of neighbours with distancees lower than the tolerance
		else if (num_of_nb > 0)
		{
			neighbours = this->get_neighbour_web(neighbours, _tol_min, _tol_max, _format);
			//remove the current atom from the neighbour list
			neighbours.erase(neighbours.begin() + (std::find(neighbours.begin(), neighbours.end(), i) - neighbours.begin()));
			//if poscar is the output format, move the neighbours which represent images to the closest image positions
			if (_format == "POSCAR")
			{
				this->move_images(i, neighbours);
			}
		}
		num_of_nb = neighbours.size();
		
		//check if there might be unreasonable "averaging" of atoms and give the user a warning
		//e.g. when three atoms get averaged to one causing the loss an atom at the border where to unit cells are combined.
		if (_format == "xyz" && num_of_nb > 2)
		{
			bad_avg = true;
		}

		else if (_format == "POSCAR" && num_of_nb > 2)
		{
			//get the average distance of the neighbours
			double tot_dist(0);
			vector3D atom_pos(this->position_matrix[i][0], this->position_matrix[i][1], this->position_matrix[i][2]);
			for (int j = 0; j < num_of_nb; j++)
			{
				int nb = neighbours[j];
				vector3D atom_pos_nb(this->position_matrix[nb][0], this->position_matrix[nb][1], this->position_matrix[nb][2]);
				tot_dist += atom_pos.calc_dist(atom_pos_nb);
			}
			tot_dist /= num_of_nb;
			if (tot_dist > 0.01)
			{
				bad_avg = true;
			}
		}
		//add the positions of the neighbours to the current atom
		for (int j = 0; j < num_of_nb; j++)
		{
			int nb = neighbours[j];
			this->position_matrix[i][0] += this->position_matrix[nb][0];
			this->position_matrix[i][1] += this->position_matrix[nb][1];
			this->position_matrix[i][2] += this->position_matrix[nb][2];
		}

		//delete the neighbours and average the position
		for (int j = 0; j < num_of_nb; j++)
		{

			int nb = neighbours[j];
			this->position_matrix.erase(this->position_matrix.begin() + (nb - j));
			this->elements_vec.erase(this->elements_vec.begin() + (nb - j));
			this->num_of_atoms--;
		}

		this->position_matrix[i][0] /= (num_of_nb + 1);
		this->position_matrix[i][1] /= (num_of_nb + 1);
		this->position_matrix[i][2] /= (num_of_nb + 1);
		i++;
	}

	if (bad_avg)
	{
		std::cout << "Warning: At least one instance found where more than two atoms are averaged to one atom. This might cause an unwanted loss of atoms. Check your final cell carefully." << std::endl;
	}
}

//finds and returns neighbours of an atom with specified index within a given range
//returns the array of the indices of the neighbours
std::vector<int> atomic_structure::find_neighbours(int atom_index, double _tol_min, double _tol_max)
{
	std::vector<int> neighbours;
	std::string el_comp_1 = this->elements_vec[atom_index];
	vector3D comp_1(this->position_matrix[atom_index][0], this->position_matrix[atom_index][1], this->position_matrix[atom_index][2]);
	for (int i = 0; i < this->num_of_atoms; i++)
	{
		std::string el_comp_2 = this->elements_vec[i];
		vector3D comp_2(this->position_matrix[i][0], this->position_matrix[i][1], this->position_matrix[i][2]);
		double d = comp_1.calc_dist(comp_2);

		if (d >= _tol_min && d <= _tol_max  && el_comp_1 == el_comp_2 && atom_index != i)
		{
			neighbours.push_back(i);
		}
	}
	return neighbours;
}

//finds neighbours of an atom with specified index within a given range considering minimum image convention
//returns the array of the indices of the neighbours
std::vector<int> atomic_structure::find_neighbours_mic(int atom_index, double _tol_min, double _tol_max)
{
	std::vector<int> neighbours;
	std::string el_comp_1 = this->elements_vec[atom_index];
	vector3D comp_1(this->position_matrix[atom_index][0], this->position_matrix[atom_index][1], this->position_matrix[atom_index][2]);
	for (int i = 0; i < this->num_of_atoms; i++)
	{
		std::string el_comp_2 = this->elements_vec[i];
		vector3D comp_2(this->position_matrix[i][0], this->position_matrix[i][1], this->position_matrix[i][2]);
		//search doubles according to minimum image convention where the images are determined by the lattice parameters.
		//if one of the images of one atom lies at a distance within the tolerance then count this atom as a neighbour
		if (el_comp_1 == el_comp_2 && atom_index != i)
		{
			bool mic_double(false);
			for (int ii = -1; ii <= 1; ii++)
			{
				for (int jj = -1; jj <= 1; jj++)
				{
					for (int kk = -1; kk <= 1; kk++)
					{
						vector3D mic_pos = comp_2 + this->lat_a * ii + this->lat_b * jj + this->lat_c * kk;
						if (comp_1.calc_dist(mic_pos) >= _tol_min && comp_1.calc_dist(mic_pos) <= _tol_max)
						{
							neighbours.push_back(i);
							mic_double = true;
							break;
						}
					}
					if (mic_double) { break; }
				}
				if (mic_double) { break; }
			}
		}
	}
	return neighbours;
}

//for a given set of neighbours, find all other neighbours recursively spanning a web of atoms where every thread has a length smaller
//than the specified range. Last input is the output format. Returns an array of all atom indices contained in the web
std::vector<int> atomic_structure::get_neighbour_web(std::vector<int> nb_web_in, double _tol_min, double _tol_max, std::string format)
{
	bool new_threads(false);
	std::vector<int> nb_web_out = nb_web_in;
	std::vector<int> new_neighbours;
	for (unsigned int i = 0; i < nb_web_in.size(); i++)
	{
		if (format == "xyz")
		{
			new_neighbours = this->find_neighbours(nb_web_in[i], _tol_min, _tol_max);
		}
		else if (format == "POSCAR")
		{
			new_neighbours = this->find_neighbours_mic(nb_web_in[i], _tol_min, _tol_max);
		}
		for (unsigned int j = 0; j < new_neighbours.size(); j++)
		{
			if (std::find(nb_web_out.begin(), nb_web_out.end(), new_neighbours[j]) == nb_web_out.end())
			{
				nb_web_out.push_back(new_neighbours[j]);
				new_threads = true;
			}
		}
	}
	if (new_threads == true)
	{
		nb_web_out = get_neighbour_web(nb_web_out, _tol_min, _tol_max, format);
	}
	std::sort(nb_web_out.begin(), nb_web_out.end());
	return nb_web_out;
}

//move the found neighbours of an atom to the image position that is closest to this atom
//the neighbour array must not include the atom of the specified index
void atomic_structure::move_images(int _atom_index, const std::vector<int>& _neighbours)
{
	vector3D comp_1(this->position_matrix[_atom_index][0], this->position_matrix[_atom_index][1], this->position_matrix[_atom_index][2]);
	for (unsigned int i = 0; i < _neighbours.size(); i++)
	{
		int nb = _neighbours[i];
		vector3D comp_2(this->position_matrix[nb][0], this->position_matrix[nb][1], this->position_matrix[nb][2]);
		double min_dist = comp_1.calc_dist(comp_2);

		for (int ii = -1; ii <= 1; ii++)
		{
			for (int jj = -1; jj <= 1; jj++)
			{
				for (int kk = -1; kk <= 1; kk++)
				{
					vector3D mic_pos = comp_2 + this->lat_a * ii + this->lat_b * jj + this->lat_c * kk;
					double dist = comp_1.calc_dist(mic_pos);
					if (dist < min_dist)
					{
						min_dist = dist;
						this->position_matrix[nb][0] = mic_pos.x;
						this->position_matrix[nb][1] = mic_pos.y;
						this->position_matrix[nb][2] = mic_pos.z;
					}
				}
			}
		}
	}
}


//sort the atoms along a specific direction and speficy if it should be sorted in blocked format
void atomic_structure::sort(const vector3D& _dir, bool _bl)
{
	std::vector<std::string> elements_vec_sorted;
	std::vector<std::vector<double>> position_matrix_sorted;
	int sorted_atoms(0);
	for (int i = 0; i < this->num_of_atoms; i++)
	{
		vector3D atom_pos(this->position_matrix[i][0], this->position_matrix[i][1], this->position_matrix[i][2]);
		if (i == 0)
		{
			elements_vec_sorted.push_back(this->elements_vec[0]);
			position_matrix_sorted.push_back(this->position_matrix[0]);
			sorted_atoms++;
		}
		bool sorted(false);
		for (int j = 0; j < sorted_atoms; j++)
		{
			vector3D atom_pos_sorted(position_matrix_sorted[j][0], position_matrix_sorted[j][1], position_matrix_sorted[j][2]);
			if (atom_pos_sorted * _dir >= atom_pos * _dir && i != 0)
			{
				position_matrix_sorted.insert(position_matrix_sorted.begin() + j, this->position_matrix[i]);
				elements_vec_sorted.insert(elements_vec_sorted.begin() + j, this->elements_vec[i]);
				sorted_atoms++;
				sorted = true;
				break;
			}
		}
		if (!sorted && i != 0)
		{
			position_matrix_sorted.push_back(this->position_matrix[i]);
			elements_vec_sorted.push_back(this->elements_vec[i]);
			sorted = true;
			sorted_atoms++;
		}
	}

	if (_bl)
	{
		//rearrange the atoms into blocks of the same atom
		std::vector<std::string> elements_vec_rear;
		std::vector<std::vector<double>> position_matrix_rear;
		for (unsigned int i = 0; i < this->elements.size(); i++)
		{
			std::string el = this->elements[i];
			for (int j = 0; j < this->num_of_atoms; j++)
			{
				if (elements_vec_sorted[j] == el)
				{
					elements_vec_rear.push_back(el);
					position_matrix_rear.push_back(position_matrix_sorted[j]);
				}
			}
		}
		this->elements_vec = elements_vec_rear;
		this->position_matrix = position_matrix_rear;
	}
	else
	{
		this->elements_vec = elements_vec_sorted;
		this->position_matrix = position_matrix_sorted;
	}
}

//set the lattice vectors
void atomic_structure::set_lat(const vector3D& _vec1, const vector3D& _vec2, const vector3D& _vec3)
{
	this->lat_a = _vec1;
	this->lat_b = _vec2;
	this->lat_c = _vec3;
}

//automatic determination of the unit cell 
//this is very crude and only possible if:
//	1. the structure is at least orthorhombic
//	2. the atomic structure contains only the atoms in the unit cell (i.e. no cuts necessary)
//	3. there is at least one atom at the 0 coordinate in every direction
void atomic_structure::find_uc()
{
	//get the maximum and minium atom cartesian coordinates in all directions
	std::vector<double> min_max_vec = this->min_max();
	double min_x = min_max_vec[0];
	double max_x = min_max_vec[1];
	double min_y = min_max_vec[2];
	double max_y = min_max_vec[3];
	double min_z = min_max_vec[4];
	double max_z = min_max_vec[5];

	//get the lattice parameters from the maximum and minium positions
	double a = max_x - min_x;
	double b = max_y - min_y;
	double c = max_z - min_z;

	//construct the unit cell using these lattice parameters. If for some reason not all atoms are in the unit cell increase
	//the lattice parameters
	int atoms_before = this->num_of_atoms;
	int atoms_after(-1);
	//arbitrary tolerance in angstrom
	double tol(0.0000001);
	atomic_structure temp = *this;
	while (atoms_before != atoms_after)	//this should be fulfilled after the first cut
	{
		*this = temp;
		vector3D lat_a(1, 0, 0, a), lat_b(0, 1, 0, b), lat_c(0, 0, 1, c);
		this->set_lat(lat_a, lat_b, lat_c);
		this->cut(-min_x, -min_y, -min_z);
		atoms_after = this->num_of_atoms;
		a += tol;
		b += tol;
		c += tol;
	}
}

//find and return the minimum and maximum cartesian coordinates of the atoms
std::vector<double> atomic_structure::min_max()
{
	std::vector<double> min_max_vec(6, 0);
	double x_min = this->position_matrix[0][0];
	double x_max = x_min;
	double y_min = this->position_matrix[0][1];
	double y_max = y_min;
	double z_min = this->position_matrix[0][2];
	double z_max = z_min;
	for (int i = 1; i < this->num_of_atoms; i++)
	{
		if (this->position_matrix[i][0] < x_min)
		{
			x_min = this->position_matrix[i][0];
		}
		else if (this->position_matrix[i][0] > x_max)
		{
			x_max = this->position_matrix[i][0];
		}
		if (this->position_matrix[i][1] < y_min)
		{
			y_min = this->position_matrix[i][1];
		}
		else if (this->position_matrix[i][1] > y_max)
		{
			y_max = this->position_matrix[i][1];
		}
		if (this->position_matrix[i][2] < z_min)
		{
			z_min = this->position_matrix[i][2];
		}
		else if (this->position_matrix[i][2] > z_max)
		{
			z_max = this->position_matrix[i][2];
		}
	}
	min_max_vec[0] = x_min;
	min_max_vec[1] = x_max;
	min_max_vec[2] = y_min;
	min_max_vec[3] = y_max;
	min_max_vec[4] = z_min;
	min_max_vec[5] = z_max;
	return min_max_vec;
}

//generate a wall in the specified direction using the current and a second atomic structures with domains of the specifieds sizes
void atomic_structure::generate_wall(const atomic_structure& _domain_2, std::string _domain_size_1, std::string _domain_size_2,
	int _wall_dir, double _tol_min, double _tol_max)
{
	//generate the first domain
	std::stringstream ss(_domain_size_1);
	int domain_size_1_x(0), domain_size_1_y(0), domain_size_1_z(0);
	ss >> domain_size_1_x >> domain_size_1_y >> domain_size_1_z;
	atomic_structure uc(*this);
	std::vector<double> T_inv = this->calc_T_inv();
	uc.complete_uc(T_inv, 0, _tol_min, _tol_max);
	atomic_structure domain_1 = *this;
	domain_1.complete_uc(T_inv, 0, _tol_min, _tol_max);

	for (int i = 0; i < domain_size_1_x; i++)
	{
		for (int j = 0; j < domain_size_1_y; j++)
		{
			for (int k = 0; k < domain_size_1_z; k++)
			{
				atomic_structure new_cell(uc);
				double shift_x = i*this->lat_a.x + j*this->lat_b.x + k*this->lat_c.x;
				double shift_y = i*this->lat_a.y + j*this->lat_b.y + k*this->lat_c.y;
				double shift_z = i*this->lat_a.z + j*this->lat_b.z + k*this->lat_c.z;
				new_cell.shift(shift_x, shift_y, shift_z);
				domain_1 = domain_1.merge(new_cell);
			}
		}
	}
	//calculate the shift for the second domain
	double domain_shift_x = domain_size_1_x * (int)(_wall_dir == 1) * this->lat_a.x +
		domain_size_1_y * (int)(_wall_dir == 2) * this->lat_b.x +
		domain_size_1_z * (int)(_wall_dir == 3) * this->lat_c.x;
	double domain_shift_y = domain_size_1_x * (int)(_wall_dir == 1) * this->lat_a.y +
		domain_size_1_y * (int)(_wall_dir == 2) * this->lat_b.y +
		domain_size_1_z * (int)(_wall_dir == 3) * this->lat_c.y;
	double domain_shift_z = domain_size_1_x * (int)(_wall_dir == 1) * this->lat_a.z +
		domain_size_1_y * (int)(_wall_dir == 2) * this->lat_b.z +
		domain_size_1_z * (int)(_wall_dir == 3) * this->lat_c.z;

	//generate the second domain
	ss.str(_domain_size_2);
	ss.clear();
	int domain_size_2_x(0), domain_size_2_y(0), domain_size_2_z(0);
	ss >> domain_size_2_x >> domain_size_2_y >> domain_size_2_z;
	uc = _domain_2;
	uc.complete_uc(T_inv, 0, _tol_min, _tol_max);
	atomic_structure domain_2 = _domain_2;
	domain_2.complete_uc(T_inv, 0, _tol_min, _tol_max);
	for (int i = 0; i < domain_size_2_x; i++)
	{
		for (int j = 0; j < domain_size_2_y; j++)
		{
			for (int k = 0; k < domain_size_2_z; k++)
			{
				atomic_structure new_cell(uc);
				double shift_x = i*this->lat_a.x + j*this->lat_b.x + k*this->lat_c.x;
				double shift_y = i*this->lat_a.y + j*this->lat_b.y + k*this->lat_c.y;
				double shift_z = i*this->lat_a.z + j*this->lat_b.z + k*this->lat_c.z;
				new_cell.shift(shift_x, shift_y, shift_z);
				domain_2 = domain_2.merge(new_cell);
			}
		}
	}

	//shift the second domain in wall direction
	domain_2.shift(domain_shift_x, domain_shift_y, domain_shift_z);

	//calculate the wall lattice parameters
	double a_1 = this->lat_a.length;
	double b_1 = this->lat_b.length;
	double c_1 = this->lat_c.length;
	double a_2 = _domain_2.lat_a.length;
	double b_2 = _domain_2.lat_b.length;
	double c_2 = _domain_2.lat_c.length;
	vector3D a_wall, b_wall, c_wall;
	switch (_wall_dir)
	{
	case 1:
	{
		a_wall = this->lat_a * domain_size_1_x + (vector3D)_domain_2.lat_a * domain_size_2_x;
		if (domain_size_1_y > domain_size_2_y || (domain_size_1_y == domain_size_2_y && b_1 >= b_2))
		{
			b_wall = this->lat_b * domain_size_1_y;
		}
		else
		{
			b_wall = (vector3D)_domain_2.lat_b * domain_size_2_y;
		}
		if (domain_size_1_z > domain_size_2_z || (domain_size_1_z == domain_size_2_z && c_1 >= c_2))
		{
			c_wall = this->lat_c * domain_size_1_z;
		}
		else
		{
			c_wall = (vector3D)_domain_2.lat_c * domain_size_2_z;
		}
		break;
	}
	case 2:
	{
		b_wall = this->lat_b * domain_size_1_y + (vector3D)_domain_2.lat_b * domain_size_2_y;
		if (domain_size_1_x > domain_size_2_x || (domain_size_1_x == domain_size_2_x && a_1 >= a_2))
		{
			a_wall = this->lat_a * domain_size_1_x;
		}
		else
		{
			a_wall = (vector3D)_domain_2.lat_a * domain_size_2_x;
		}
		if (domain_size_1_z > domain_size_2_z || (domain_size_1_z == domain_size_2_z && c_1 >= c_2))
		{
			c_wall = this->lat_c * domain_size_1_z;
		}
		else
		{
			c_wall = (vector3D)_domain_2.lat_c * domain_size_2_z;
		}
		break;
	}
	case 3:
	{
		c_wall = this->lat_c * domain_size_1_z + (vector3D)_domain_2.lat_c * domain_size_2_z;
		if (domain_size_1_x > domain_size_2_x || (domain_size_1_x == domain_size_2_x && a_1 >= a_2))
		{
			a_wall = this->lat_a * domain_size_1_x;
		}
		else
		{
			a_wall = (vector3D)_domain_2.lat_a * domain_size_2_x;
		}
		if (domain_size_1_y > domain_size_2_y || (domain_size_1_y == domain_size_2_y && b_1 >= b_2))
		{
			b_wall = this->lat_b * domain_size_1_y;
		}
		else
		{
			b_wall = (vector3D)_domain_2.lat_b * domain_size_2_y;
		}
		break;
	}
	default:
		break;
	}

	//merge the domains
	*this = domain_1.merge(domain_2);

	//set the new lattice parameters
	this->lat_a = a_wall;
	this->lat_b = b_wall;
	this->lat_c = c_wall;
}

//recursively find and add all missing atoms on the unit cell boundaries that are lost when reading a poscar file.
//inputs are the base change matrix to get the relative coordinates and the index of the atom for which images on the boundaries are looked for
void atomic_structure::complete_uc(const std::vector<double>& _T_inv, int _index, double _tol_min, double _tol_max)
{
	//tolerances of what is considered to still lie on the unit cell boundaries
	double x_tol = 1 / this->lat_a.length;
	double y_tol = 1 / this->lat_b.length;
	double z_tol = 1 / this->lat_c.length;


	if (_index < this->num_of_atoms)
	{
		std::string el = this->elements_vec[_index];
		vector3D atom_pos(this->position_matrix[_index][0], this->position_matrix[_index][1], this->position_matrix[_index][2]);
		vector3D image;
		std::vector<vector3D> image_vec(0);
		double x_direct = _T_inv[0] * atom_pos.x + _T_inv[1] * atom_pos.y + _T_inv[2] * atom_pos.z;
		double y_direct = _T_inv[3] * atom_pos.x + _T_inv[4] * atom_pos.y + _T_inv[5] * atom_pos.z;
		double z_direct = _T_inv[6] * atom_pos.x + _T_inv[7] * atom_pos.y + _T_inv[8] * atom_pos.z;

		x_direct = abs(x_direct);
		y_direct = abs(y_direct);
		z_direct = abs(z_direct);
		//find all images of the current atom
		if (x_direct <= x_tol)
		{
			image = atom_pos + this->lat_a;
			image_vec.push_back(image);
		}
		if (abs(1 - x_direct) <= x_tol)
		{
			image = atom_pos - this->lat_a;
			image_vec.push_back(image);
		}
		if (y_direct <= y_tol)
		{
			image = atom_pos + this->lat_b;
			image_vec.push_back(image);
		}
		if (abs(1 - y_direct) <= y_tol)
		{
			image = atom_pos - this->lat_b;
			image_vec.push_back(image);
		}
		if (z_direct <= z_tol)
		{
			image = atom_pos + this->lat_c;
			image_vec.push_back(image);
		}
		if (abs(1 - z_direct) <= z_tol)
		{
			image = atom_pos - this->lat_c;
			image_vec.push_back(image);
		}
		//check if at least one image was found
		if (image_vec.size() > 0)
		{
			//check all images
			for (unsigned int i = 0; i < image_vec.size(); i++)
			{
				//check if there is already an atom at image position
				bool occupied(false);
				vector3D image_temp(image_vec[i].x, image_vec[i].y, image_vec[i].z);
				for (int j = 0; j < this->num_of_atoms; j++)
				{
					vector3D atom_pos_temp(this->position_matrix[j][0], this->position_matrix[j][1], this->position_matrix[j][2]);
					double dist = image_temp.calc_dist(atom_pos_temp);
					if (dist >= _tol_min && dist <= _tol_max && _tol_min != 0 && _tol_max != 0)
					{
						occupied = true;
						break;
					}
				}
				//if there isn't, add the image
				if (!occupied)
				{
					this->elements_vec.push_back(el);
					this->position_matrix.push_back({ image_temp.x, image_temp.y, image_temp.z });
					this->num_of_atoms++;
				}
			}
		}
		_index++;
		this->complete_uc(_T_inv, _index, _tol_min, _tol_max);
	}
	else
	{
		return;
	}
}

//caluclates and returns the base change matrix to get from the cartesian to relative coordinates as determined by the lattice vectors
std::vector<double> atomic_structure::calc_T_inv()
{
	std::vector<double> T(9, 0);
	std::vector<double> T_inv(9, 0);
	T[0] = this->lat_a.x;
	T[1] = this->lat_b.x;
	T[2] = this->lat_c.x;
	T[3] = this->lat_a.y;
	T[4] = this->lat_b.y;
	T[5] = this->lat_c.y;
	T[6] = this->lat_a.z;
	T[7] = this->lat_b.z;
	T[8] = this->lat_c.z;

	double det_T = T[0] * (T[4] * T[8] - T[5] * T[7])
		- T[1] * (T[3] * T[8] - T[5] * T[6])
		+ T[2] * (T[3] * T[7] - T[4] * T[6]);

	T_inv[0] = T[4] * T[8] - T[5] * T[7];
	T_inv[1] = T[2] * T[7] - T[1] * T[8];
	T_inv[2] = T[1] * T[5] - T[2] * T[4];
	T_inv[3] = T[5] * T[6] - T[3] * T[8];
	T_inv[4] = T[0] * T[8] - T[2] * T[6];
	T_inv[5] = T[2] * T[3] - T[0] * T[5];
	T_inv[6] = T[3] * T[7] - T[4] * T[6];
	T_inv[7] = T[1] * T[6] - T[0] * T[7];
	T_inv[8] = T[0] * T[4] - T[1] * T[3];

	std::transform(T_inv.begin(), T_inv.end(), T_inv.begin(), std::bind1st(std::multiplies<double>(), 1 / det_T));

	return T_inv;
}

//_________________________________________________________________________________________________________//

int main()
{

	//check inputs
	std::ifstream ifs1(input_file_1), ifs2(input_file_2);
	std::string input_format_1 = input_file_1.substr(input_file_1.find_last_of(".") + 1);
	std::string input_format_2 = input_file_2.substr(input_file_2.find_last_of(".") + 1);
	if (!ifs1.good())
	{
		std::cerr << "Error: 1st input file does not exist." << std::endl;
		keep_window_open(open_window);
		return 1;
	}

	if (!ifs2.good() && !input_file_2.empty())
	{
		std::cerr << "Error: 2nd input file does not exist." << std::endl;
		keep_window_open(open_window);
		return 1;
	}

	if (output_format != "POSCAR" && output_format != "xyz")
	{
		std::cerr << "Error: Output format not supported. Choose either \"POSCAR\" or \"xyz\" format." << std::endl;
		keep_window_open(open_window);
		return 1;
	}
	if ((input_format_1 != "xyz" && input_format_1 != "vasp")
		|| input_format_1 == input_file_1)
	{
		std::cout << "Warning: file extension of 1st input file not recognized (input files need to be in .xyz or .vasp format)" << std::endl;
		std::cout << "Trying to read as .vasp file." << std::endl << std::endl;
		input_format_1 = "POSCAR";
	}
	if ((input_format_2 != "xyz" && input_format_2 != "vasp")
		|| input_format_2 == input_file_2)
	{
		std::cout << "Warning: file extension of 2nd input file not recognized (input files need to be in .xyz or .vasp format)" << std::endl;
		std::cout << "Trying to read as .vasp file." << std::endl << std::endl;
		input_format_2 = "POSCAR";
	}

	if (!auto_uc && !wall_uc && (a == 0 || b == 0 || c == 0))
	{
		if (output_format == "xyz")
		{
			std::cout << "Warning: At least one lattice parameter is 0." << std::endl;
		}
		else if (output_format == "POSCAR" && (input_format_1 != "POSCAR" || input_format_2 != "POSCAR"))
		{
			std::cerr << "Error: No proper unit cell defined." << std::endl;
			keep_window_open(open_window);
			return 1;
		}
	}

	if (block == false && output_format == "POSCAR")
	{
		std::cerr << "Error: BLOCK needs to be true for POSCAR output" << std::endl;
		keep_window_open(open_window);
		return 1;
	}

	if (wall_dir < 1 || wall_dir > 3)
	{
		std::cerr << "Error: WALL_DIR needs to be 1, 2 or 3" << std::endl;
		keep_window_open(open_window);
		return 1;
	}
	std::stringstream check_domain_size(domain_size_1);
	int num(0);
	while (check_domain_size >> num || !check_domain_size.eof())
	{
		if (check_domain_size.fail())
		{
			std::cerr << "Error: syntax for 1st domain wall size wrong" << std::endl;
			keep_window_open(open_window);
			return 1;
		}
	}
	check_domain_size.str(domain_size_2);
	check_domain_size.clear();
	while (check_domain_size >> num || !check_domain_size.eof())
	{
		if (check_domain_size.fail())
		{
			std::cerr << "Error: syntax for 2nd domain wall size wrong" << std::endl;
			keep_window_open(open_window);
			return 1;
		}
	}

	if ((input_format_1 == "xyz" || input_format_2 == "xyz") && wall_uc)
	{
		std::cerr << "Error: to generate domain walls, input files need to be in POSCAR format" << std::endl;
		keep_window_open(open_window);
		return 1;
	}

	//check file integrity
	try
	{
		atomic_structure file_1(input_file_1, input_format_1);
	}
	catch (const std::string& exc)
	{
		std::cerr << exc << " Please check the 1st input file." << std::endl;
		keep_window_open(open_window);
		return 1;
	}

	try
	{
		atomic_structure file_2(input_file_2, input_format_2);
	}
	catch (const std::string& exc)
	{
		if (ifs2.good() && !input_file_2.empty())
		{
			std::cerr << exc << " Please check the 2nd input file." << std::endl;
			keep_window_open(open_window);
			return 1;
		}
	}

	atomic_structure file_1(input_file_1, input_format_1);
	atomic_structure file_2;
	if (ifs2.good() && !input_file_2.empty())
	{
		atomic_structure file_2_temp(input_file_2, input_format_2);
		file_2 = file_2_temp;
	}

	vector3D rot_a(rot_a_x, rot_a_y, rot_a_z, 1);
	vector3D rot_b(rot_b_x, rot_b_y, rot_b_z, 1);
	vector3D rot_c(rot_c_x, rot_c_y, rot_c_z, 1);
	if (alpha1 != 0 || beta1 != 0 || gamma1 != 0)
	{
		std::cout << std::endl << "rotating 1st cell by: " << std::endl;
		std::cout << "\t" << alpha1 << " degrees around [" << rot_a_x << "," << rot_a_y << "," << rot_a_z << "]" << std::endl;
		std::cout << "\t" << beta1 << " degrees around [" << rot_b_x << "," << rot_b_y << "," << rot_b_z << "]" << std::endl;
		std::cout << "\t" << gamma1 << " degrees around [" << rot_c_x << "," << rot_c_y << "," << rot_c_z << "]" << std::endl;
		file_1.rotate(alpha1, beta1, gamma1, rot_a, rot_b, rot_c);
	}
	if ((alpha2 != 0 || beta2 != 0 || gamma2 != 0) && !input_file_2.empty())
	{
		std::cout << std::endl << "rotating 2nd cell by: " << std::endl;
		std::cout << "\t" << alpha2 << " degrees around [" << rot_a_x << "," << rot_a_y << "," << rot_a_z << "]" << std::endl;
		std::cout << "\t" << beta2 << " degrees around [" << rot_b_x << "," << rot_b_y << "," << rot_b_z << "]" << std::endl;
		std::cout << "\t" << gamma2 << " degrees around [" << rot_c_x << "," << rot_c_y << "," << rot_c_z << "]" << std::endl;
		file_2.rotate(alpha2, beta2, gamma2, rot_a, rot_b, rot_c);
	}

	std::vector<double> min_max_1 = file_1.min_max();
	std::cout << std::endl << "minimum and maximum positions of 1st cell before merging: " << std::endl;
	std::cout << '\t' << "X: [" << min_max_1[0] << "," << min_max_1[1] << "]" << std::endl;
	std::cout << '\t' << "Y: [" << min_max_1[2] << "," << min_max_1[3] << "]" << std::endl;
	std::cout << '\t' << "Z: [" << min_max_1[4] << "," << min_max_1[5] << "]" << std::endl;

	if (!input_file_2.empty())
	{
		std::vector<double> min_max_2 = file_2.min_max();
		std::cout << std::endl << "minimum and maximum positions of 2nd cell before merging: " << std::endl;
		std::cout << '\t' << "X: [" << min_max_2[0] << "," << min_max_2[1] << "]" << std::endl;
		std::cout << '\t' << "Y: [" << min_max_2[2] << "," << min_max_2[3] << "]" << std::endl;
		std::cout << '\t' << "Z: [" << min_max_2[4] << "," << min_max_2[5] << "]" << std::endl;
	}

	if (shift_x_1 != 0 || shift_y_1 != 0 || shift_z_1 != 0)
	{
		std::cout << std::endl << "shifting 1st cell by [" << shift_x_1 << "," << shift_y_1 << "," << shift_z_1 << "]" << std::endl;
		file_1.shift(shift_x_1, shift_y_1, shift_z_1);
	}

	if ((shift_x_2 != 0 || shift_y_2 != 0 || shift_z_2 != 0) && !input_file_2.empty())
	{
		std::cout << std::endl << "shifting 2nd cell by [" << shift_x_2 << "," << shift_y_2 << "," << shift_z_2 << "]" << std::endl;
		file_2.shift(shift_x_2, shift_y_2, shift_z_2);
	}

	atomic_structure final_file(file_1);
	if (wall_uc)
	{
		std::cout << std::endl << "generating domain wall" << std::endl;
		final_file.generate_wall(file_2, domain_size_1, domain_size_2, wall_dir, tol_min, tol_max);
	}
	else
	{
		if (!input_file_2.empty())
		{
			std::cout << std::endl << "merging structures" << std::endl;
			final_file = final_file.merge(file_2);
			std::vector<double> min_max = final_file.min_max();
			std::cout << std::endl << "minimum and maximum positions of merged file: " << std::endl;
			std::cout << '\t' << "X: [" << min_max[0] << "," << min_max[1] << "]" << std::endl;
			std::cout << '\t' << "Y: [" << min_max[2] << "," << min_max[3] << "]" << std::endl;
			std::cout << '\t' << "Z: [" << min_max[4] << "," << min_max[5] << "]" << std::endl;
		}

		if (!auto_uc && (a != 0 || b != 0 || c != 0))
		{
			std::cout << std::endl << "cutting out atoms" << std::endl;
			vector3D lat_a(a_x, a_y, a_z, a), lat_b(b_x, b_y, b_z, b), lat_c(c_x, c_y, c_z, c);
			final_file.set_lat(lat_a, lat_b, lat_c);
			final_file.cut(o_x, o_y, o_z);
		}
		else if (auto_uc)
		{
			std::cout << std::endl << "finding unit cell" << std::endl;
			final_file.find_uc();
		}
	}

	//sorting file
	if (sort_x != 0 || sort_y != 0 || sort_z != 0)
	{
		std::cout << std::endl << "sorting atoms along [" << sort_x << "," << sort_y << "," << sort_z << "] direction" << std::endl;
		vector3D sort_dir(sort_x, sort_y, sort_z);
		final_file.sort(sort_dir, block);
	}

	//averaging doubles
	if (tol_min < tol_max)
	{
		if (output_format == "xyz")
		{
			std::cout << std::endl << "averaging all atoms within a distance of " << tol_min << " to " << tol_max << " Angstrom" << std::endl;
			final_file.average_atoms(tol_min, tol_max, output_format);
		}
		else if (output_format == "POSCAR")
		{
			std::cout << std::endl << "averaging all atoms within a distance of " << tol_min << " to " << tol_max << " Angstrom" << std::endl;
			final_file.average_atoms(tol_min, tol_max, output_format);
		}
	}

	//writing file
	if (output_format == "xyz")
	{
		std::cout << std::endl << "writing .xyz file" << std::endl;
		final_file.write_xyz(output_file);
	}
	else if (output_format == "POSCAR")
	{
		std::cout << std::endl << "writing POSCAR file" << std::endl;
		final_file.write_poscar(output_file);
	}
	keep_window_open(open_window);
}

