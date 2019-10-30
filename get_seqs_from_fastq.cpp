#include "get_seqs_from_fastq.h"
using namespace std;

int main()
{

// Open fastq file

	char* file_in = new char[255];

//	Change the string below to point to your FASTQ file

strcpy(file_in, "F:\\My Documents\\Students\\Mzwanele Ngubo\\PhD\\Sequence data\\Sample_H1\\H1_ATCACG_L001_R1_001.fastq");
cout << "Fastq file: " <<  file_in << endl;

// Read all the sequences to an array
HP_DynamicStringArray* string_array = new HP_DynamicStringArray();
HP_ReadTextFile* readfile = new HP_ReadTextFile;

int entries = readfile->ReadFastqFile(file_in, string_array);
cout << entries << " sequences retrieved from fastq file" << endl;
char** primers = new char*[8];
for(int x = 0; x < 8; x++)
{
	primers[x] = new char[18];
	memset(primers[x], 0, 18);
}

strcpy(primers[0],"ATGTCCACGAGGTCTCT");
strcpy(primers[4], "TACGCTGCAGGTCGAGG");

strcpy(primers[1], "CCTCGACCTGCAGCGTA");
strcpy(primers[5], "AGAGACCTCGTGGACAT");

strcpy(primers[2], "CGGTGTCGGTCTCGTAG");
strcpy(primers[6], "GATGAATTCGAGCTGGG");

strcpy(primers[3], "CCCAGCTCGAATTCATC");
strcpy(primers[7], "CTACGAGACCGACACCG");


int two_primers = 0;
int correct_length = 0;
int one_primer = 0;
int incorrect_length = 0;
int no_primers = 0;

// See if we can match a primer

int primer_length = 17;
int insert_length = 20;
char* pointer_1;
char* pointer_2;
int primer_1;
int primer_2;
HP_DynamicStringArray* barcodes = new HP_DynamicStringArray();
HPDynamicIntArray* up_down = new HPDynamicIntArray();
char* temp = new char[21];
for(int x = 0; x < string_array->GetNumberOfStrings(); x++)
	{
	memset(temp, 0, 21);
	int index = 0;
	while(index < 4 && !strstr(string_array->GetStringPointer(x), primers[index]))
		index++;
	if(index < 4) // we  found a match
		{
		pointer_1 = strstr(string_array->GetStringPointer(x), primers[index]);
		// See if we can find the match for the primer pair
		pointer_2 = strstr(string_array->GetStringPointer(x), primers[index+4]);
		primer_1 = pointer_1 - string_array->GetStringPointer(x);
		primer_2 = pointer_2 - string_array->GetStringPointer(x);
		if(pointer_2 == 0) // could not find a match
			{
				one_primer++;
			}
		else
			{
			if(primer_2 > primer_1) // the order of primers is correct
				{
				two_primers++;
				if((primer_2 - primer_1 - primer_length) == insert_length)
					{
						correct_length++;
						strncpy(temp, (pointer_1+primer_length), insert_length);
						if(index < 2)
							up_down->AddInt(1);
						else
							up_down->AddInt(0);
						if((index == 1) || (index == 3))
							barcodes->AddString(readfile->ReverseComplement(temp));
						else
							barcodes->AddString(temp);


					}
				else
					incorrect_length++;
				}
			}

		}
	else
		{
		no_primers++;
		}
	}

cout << "two_primers: " << two_primers << endl;
cout << "correct_length: " << correct_length << endl;
cout << "one_primer: " << one_primer << endl;
cout << "incorrect_length: " << incorrect_length << endl;
cout << "no_primers: " << no_primers << endl;

// Save the strings from the array

char* outfilename = new char[1024];

//	Change the string below to a filename and path of your choice

strcpy(outfilename, "F:\\My Documents\\Students\\Mzwanele Ngubo\\PhD\\Sequence data\\Sample_H1\\H1_ATCACG_L001_R1_001_barcodes.txt");
fstream outfile;
outfile.open(outfilename, fstream::out | fstream::binary);
char* temp_2 = new char[30];
for(int x = 0; x < barcodes->GetNumberOfStrings(); x++)
{
	memset(temp_2, 0, 30);
	sprintf(temp_2, "%s%c%d", barcodes->GetStringPointer(x), '\t', up_down->GetEntry(x));
	outfile << temp_2 << endl;
}
outfile.close();

cout << "done" << endl;

return 0;
}
