/* Calculate ODF C++ Version
By Edmund Wong
Current variables used: imvol, qdata, odf_faces, odf_vertices
*/

//KNOWN ISSUES: nifti_image_read_bricks only reads images with 128x128x128x515 volumes properly.
// This is possibly a bug with nifti_image_read_bricks
// For example, it reads the seal image (64x64x32x515)  correctly for the first 6291456 elements, then after that, the values are wrong
// Needs user input for path to find the nifti_file.
// Needs user input for output name

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "nifti1_io.h"
#include <math.h>
#include "fftw3.h"

int main()
{
    using namespace std;
    
    //Timer start
    time_t rawtime;
    struct tm * timeinfo;
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    
    ifstream in, in2, in3;
    string qdatapath, qdirfile, qfacefile, qverticefile;

    //Prepare user input
    while (true)
    {
        cout << "Please specify the path of the necessary ODF processing files." << endl;
        cout << "Path needs to contain q.csv, odf_faces.csv, and odf_vertices.csv." << endl;
        cin >> qdatapath;
        qdirfile = qdatapath;
        qfacefile = qdatapath;
        qverticefile = qdatapath;
        
        qdirfile.append("q.csv");
        qfacefile.append("odf_faces.csv");
        qverticefile.append("odf_vertices.csv");
        
        char *qfilename1=new char[qdirfile.size()+1];
        char *qfilename2=new char[qfacefile.size()+1];
        char *qfilename3=new char[qverticefile.size()+1];
        qfilename1[qdirfile.size()]=0;
        qfilename2[qfacefile.size()]=0;
        qfilename3[qverticefile.size()]=0;
        memcpy(qfilename1,qdirfile.c_str(),qdirfile.size());
        memcpy(qfilename2,qfacefile.c_str(),qfacefile.size());
        memcpy(qfilename3,qverticefile.c_str(),qverticefile.size());
        
        in.open(qfilename1);
        in2.open(qfilename2);
        in3.open(qfilename3);
        if ((!in.good()) || (!in2.good()) || (!in3.good()))
        {
            in.close();
            in2.close();
            in3.close();
            cout << "**ERROR: The path is not valid or the path does not contain all of the necessary odf processing files." << endl << endl;
        }
        else
        {
            cout << "Loading " << qdirfile << endl;
            cout << "Loading " << qfacefile << endl;
            cout << "Loading " << qverticefile << endl;
            break;
        }
    }
    
    cout << "START TIME/DATE: " << asctime (timeinfo) << endl;
    
    string line, field;
    int qdata[515][3];
    string s;
    char *a;
    
    // READ IMAGE
    nifti_image *img;
    nifti_brick_list NB_orig;
    signed int *imvol; //change to float or int depending on dsi raw image
    
    img = nifti_image_read_bricks("/scratch/my_dsi_file.nii", 0, NULL, &NB_orig);
    imvol=(signed int*)NB_orig.bricks[0]; //change to float or short depending on dsi raw image

    int xdim = img->dim[1];
    int ydim = img->dim[2];
    int zdim = img->dim[3];
    int tdim = img->dim[4];
    
    vector< vector<string> > array;  // the 2D array
    vector<string> v;                // array of values for one line only
    vector< vector<string> > array2;  // the 2D array
    vector<string> v2;
    vector< vector<string> > array3;  // the 2D array
    vector<string> v3;

    //Read in the data from file one (q.csv) and arrange into an array
    
    while ( getline(in,line) )
    {
        v.clear();
        stringstream ss(line);        
        while (getline(ss,field,','))
            v.push_back(field);        
        array.push_back(v);
    }
    
    
    for (size_t i=0; i<array.size(); ++i)
        for (size_t j=0; j<array[i].size(); ++j)
        {
            s = array[i][j];
            a=new char[s.size()+1];
            a[s.size()]=0;
            memcpy(a,s.c_str(),s.size());
            qdata[i][j]=atoi(a);
        }
    
    //Read in the data from file one (ODF FACES) and arrange into an array
    
    while ( getline(in2,line) )
    {
        v2.clear();
        stringstream ss(line);
        while (getline(ss,field,','))
            v2.push_back(field);
        array2.push_back(v2);
    }
    
    int no_odf_faces = (int)array2[0].size();
    int odf_faces[3][no_odf_faces];
    
    for (size_t i=0; i<array2.size(); ++i)
        for (size_t j=0; j<array2[i].size(); ++j)
        {
            s = array2[i][j];
            a=new char[s.size()+1];
            a[s.size()]=0;
            memcpy(a,s.c_str(),s.size());
            odf_faces[i][j]=atoi(a);
        }
    
    //Read in the data from file one (ODF VERTICES) and arrange into an array
    
    while ( getline(in3,line) )
    {
        v3.clear();
        stringstream ss(line);        
        while (getline(ss,field,','))
            v3.push_back(field);        
        array3.push_back(v3);
    }
    
    int no_odf_vertices = (int)array3[0].size();
    int half_no_odf_vertices = no_odf_vertices/2;
    double odf_vertices[3][no_odf_vertices];
    
    for (size_t i=0; i<array3.size(); ++i)
        for (size_t j=0; j<array3[i].size(); ++j)
        {
            s = array3[i][j];
            a=new char[s.size()+1];
            a[s.size()]=0;
            memcpy(a,s.c_str(),s.size());
            odf_vertices[i][j]=atof(a);
        }
    
    // Declaring variables
    double hanning_filter[515];
    double value[tdim];    
    double sq[16][16][16];
    double sqnew[16][16][16];    
    double xi[20], yi[20], zi[20], r[20], v_xyz[20], x_val, y_val, z_val;
    int x_fl, x_ce, y_fl, y_ce, z_fl, z_ce;
    double origin = 9;
    int N = 16*16*16;    
    size_t ii, jj, kk;

    fftw_complex *in_fft, *out_fft;
    in_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N); // N is 16*16*16
    out_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_plan p;
    double fft_output[16][16][16];
    double Pr[16][16][16];
    double odf[half_no_odf_vertices];
    float *odf_whole = new float[xdim*ydim*zdim*half_no_odf_vertices];/*modified*/
    
    // Assign constant values to r vector
    for (size_t idx=0; idx<20; idx++)
        r[idx]=2.1+0.2*idx;
    
    // Defining the hanning filter and adjusting the qdata coordinates
    for (size_t idx=0; idx<515; idx++)
    {
        hanning_filter[idx] = 0.5*(1+cos(2*PI*sqrt(qdata[idx][0]*qdata[idx][0]+qdata[idx][1]*qdata[idx][1]+qdata[idx][2]*qdata[idx][2])/16));
        qdata[idx][0] = qdata[idx][0]+8; //9-1
        qdata[idx][1] = qdata[idx][1]+8;
        qdata[idx][2] = qdata[idx][2]+8;
    }
    
    for (size_t x_idx=0; x_idx<xdim; x_idx++)
    {
        cout << "Progress: " << x_idx << " of 127" << endl;
        for (size_t y_idx=0; y_idx<ydim; y_idx++)
        {
            for (size_t z_idx=0; z_idx<zdim; z_idx++)
            {
                // Initialize a 16*16*16 space with zeros
                for (size_t i=0; i<16; i++)
                    for (size_t j=0; j<16; j++)
                        for (size_t k=0; k<16; k++)
                            sq[i][j][k]=0;
                
                // Applying the hanning filter to the signal values, then assigning them to qspace
                for (size_t idx=0; idx<515; idx++)
                {
                    value[idx] = imvol[x_idx + y_idx*xdim + z_idx*xdim*ydim + idx*xdim*ydim*zdim] * hanning_filter[idx]; // IDENTIFY VOXEL LOCATION
                    sq[qdata[idx][0]][qdata[idx][1]][qdata[idx][2]]=value[idx];
                }
                
                // fftshift of qspace data
                for (size_t i=0; i<16; i++)
                    for (size_t j=0; j<16; j++)
                        for (size_t k=0; k<16; k++)
                        {
                            if (i>7)
                                ii=i-8;
                            else
                                ii=i+8;
                            if (j>7)
                                jj=j-8;
                            else
                                jj=j+8;
                            if (k>7)
                                kk=k-8;
                            else
                                kk=k+8;
                            sqnew[ii][jj][kk]=sq[i][j][k];
                        }
                
                // Fast Fourier Transform of qspace data
                for (size_t i=0; i<16; i++)
                    for (size_t j=0; j<16; j++)
                        for (size_t k=0; k<16; k++)
                        {
                            in_fft[k*16*16+j*16+i][0]=sqnew[i][j][k];
                        }
                p = fftw_plan_dft_3d(16,16,16, in_fft, out_fft, FFTW_FORWARD, FFTW_ESTIMATE);
                fftw_execute(p);
                
                // Assign 1-D fft output results to a 3D array
                for (size_t i=0; i<16; i++)
                    for (size_t j=0; j<16; j++)
                        for (size_t k=0; k<16; k++)
                        {
                            fft_output[i][j][k]=fabs(out_fft[k*16*16+j*16+i][0]);
                        }
                
                // fftshift again
                for (size_t i=0; i<16; i++)
                    for (size_t j=0; j<16; j++)
                        for (size_t k=0; k<16; k++)
                        {
                            if (i>7)
                                ii=i-8;
                            else
                                ii=i+8;
                            if (j>7)
                                jj=j-8;
                            else
                                jj=j+8;
                            if (k>7)
                                kk=k-8;
                            else
                                kk=k+8;
                            Pr[ii][jj][kk]=fft_output[i][j][k];
                        }
                
                // zero-lize the odf array
                for (size_t idx=0; idx<half_no_odf_vertices; idx++)
                    odf[idx]=0;
                
                for (size_t mm=0; mm<half_no_odf_vertices; mm++)
                {
                    for (int iii=0;iii<20;iii++)
                    {
                        xi[iii]=origin+r[iii]*odf_vertices[0][mm];
                        yi[iii]=origin+r[iii]*odf_vertices[1][mm];
                        zi[iii]=origin+r[iii]*odf_vertices[2][mm];
                        // Trilinear interpolation of point at (xi,yi,zi)
                        x_fl=floor(xi[iii])-1; x_ce=ceil(xi[iii])-1;
                        y_fl=floor(yi[iii])-1; y_ce=ceil(yi[iii])-1;
                        z_fl=floor(zi[iii])-1; z_ce=ceil(zi[iii])-1;
                        if (x_fl==x_ce) x_ce=x_ce+1;
                        if (y_fl==y_ce) y_ce=y_ce+1;
                        if (z_fl==z_ce) z_ce=z_ce+1;
                        x_val=xi[iii]-x_fl-1; y_val=yi[iii]-y_fl-1; z_val=zi[iii]-z_fl-1;
                        
                        // MATLAB interp3 method
                        // Somehow the results only match matlabs results if the x_fl, x_ce becomes y_fl, y_ce respectively, and vice versa
                        // It can be investigated inside the interp3 function in matlab, perhaps in lines 355 and 374. The below formula works though.
                        v_xyz[iii]=((Pr[y_fl][x_fl][z_fl]*(1-y_val)+Pr[y_ce][x_fl][z_fl]*y_val)*(1-x_val)+(Pr[y_fl][x_ce][z_fl]*(1-y_val)+Pr[y_ce][x_ce][z_fl]*y_val)*x_val)*(1-z_val)+((Pr[y_fl][x_fl][z_ce]*(1-y_val)+Pr[y_ce][x_fl][z_ce]*y_val)*(1-x_val)+(Pr[y_fl][x_ce][z_ce]*(1-y_val)+Pr[y_ce][x_ce][z_ce]*y_val)*x_val)*z_val;
                        
                        // Another method but results are slightly different but close (?), based on wikipedia and another website: http://mipav.cit.nih.gov/documentation/HTML%20Algorithms/InterplationMethods.html
                        /*
                         v_xyz2[iii]=Pr[x_fl][y_fl][z_fl]*(1-x_val)*(1-y_val)*(1-z_val)+Pr[x_ce][y_fl][z_fl]*x_val*(1-y_val)*(1-z_val)+Pr[x_fl][y_ce][z_fl]*(1-x_val)*y_val*(1-z_val)+Pr[x_fl][y_fl][z_ce]*(1-x_val)*(1-y_val)*z_val+Pr[x_ce][y_fl][z_ce]*x_val*(1-y_val)*z_val+Pr[x_fl][y_ce][z_ce]*(1-x_val)*y_val*z_val+Pr[x_ce][y_ce][z_fl]*x_val*y_val*(1-z_val)+Pr[x_ce][y_ce][z_ce]*x_val*y_val*z_val;
                         */
                        
                        odf[mm]=odf[mm]+v_xyz[iii]*pow(r[iii],2);
                    }
                    //odf_whole[x_idx+xdim*y_idx+xdim*ydim*z_idx+xdim*ydim*zdim*mm]=odf[mm];
                    odf_whole[x_idx+xdim*y_idx+xdim*ydim*z_idx+xdim*ydim*zdim*mm]=odf[mm]; //modified
                }
            }
        }
    }
    
    FILE * pFile;
    
    pFile = fopen ( "/scratch/output.nii" , "wb" );
    
    int * sizeof_hdr = new int;
    char * data_type = new char[10];
    char * db_name = new char[18];
    int * extents = new int;
    short * session_error = new short;
    char * regular = new char;
    char * dim_info = new char;
    short * dim = new short[8];
    float * intent_p1 = new float;
    float * intent_p2 = new float;
    float * intent_p3 = new float;
    short * intent_code = new short;
    short * datatype = new short;
    short * bitpix = new short;
    short * slice_start = new short;
    float * pixdim = new float[8];
    float * vox_offset = new float;
    float * scl_slope = new float;
    float * scl_inter = new float;
    short * slice_end = new short;
    char * slice_code = new char;
    char * xyzt_units = new char;
    float * cal_max = new float;
    float * cal_min = new float;
    float * slice_duration = new float;
    float * toffset = new float;
    int * glmax = new int;
    int * glmin = new int;
    char * descrip = new char[80];
    char * aux_file = new char[24];
    short * qform_code = new short;
    short * sform_code = new short;
    float * quatern_b = new float;
    float * quatern_c = new float;
    float * quatern_d = new float;
    float * qoffset_x = new float;
    float * qoffset_y = new float;
    float * qoffset_z = new float;
    float * srow_x = new float[4];
    float * srow_y = new float[4];
    float * srow_z = new float[4];
    char * intent_name = new char[16];
    char * magic = new char[4];
    * sizeof_hdr = 348;
    * regular = 114;
    // Dimensions
    dim[0] = 4; dim[1]=xdim; dim[2]=ydim; dim[3]=zdim; dim[4]=half_no_odf_vertices; dim[5]=1; dim[6]=1; dim[7]=1;
    * datatype = 16;
    * bitpix = 64;
    pixdim[0] = 0; pixdim[1]=0.4688; pixdim[2]=0.4688; pixdim[3]=0.6; pixdim[4]=45000;
    * vox_offset = 352;
    * xyzt_units = 10;
    magic[0]=110; magic[1]=43; magic[2]=49;
    
    fwrite (sizeof_hdr, sizeof(int), 1, pFile);
    fwrite (data_type, sizeof(char), 10, pFile);
    fwrite (db_name, sizeof(char), 18, pFile);
    fwrite (extents, sizeof(int), 1, pFile);
    fwrite (session_error, sizeof(short), 1, pFile);
    fwrite (regular, sizeof(char), 1, pFile);
    fwrite (dim_info, sizeof(char), 1, pFile);
    fwrite (dim, sizeof(short), 8, pFile);
    fwrite (intent_p1, sizeof(float), 1, pFile);
    fwrite (intent_p2, sizeof(float), 1, pFile);
    fwrite (intent_p3, sizeof(float), 1, pFile);
    fwrite (intent_code, sizeof(short), 1, pFile);
    fwrite (datatype, sizeof(short),1, pFile);
    fwrite (bitpix, sizeof(short), 1, pFile);
    fwrite (slice_start, sizeof(short), 1, pFile);
    fwrite (pixdim, sizeof(float), 8, pFile);
    fwrite (vox_offset, sizeof(float), 1, pFile);
    fwrite (scl_slope, sizeof(float), 1, pFile);
    fwrite (scl_inter, sizeof(float), 1, pFile);
    fwrite (slice_end, sizeof(short), 1, pFile);
    fwrite (slice_code, sizeof(char), 1, pFile);
    fwrite (xyzt_units, sizeof(char), 1, pFile);
    fwrite (cal_max, sizeof(float), 1, pFile);
    fwrite (cal_min, sizeof(float), 1, pFile);
    fwrite (slice_duration, sizeof(float), 1, pFile);
    fwrite (toffset, sizeof(float), 1, pFile);
    fwrite (glmax, sizeof(int), 1, pFile);
    fwrite (glmin, sizeof(int), 1, pFile);
    fwrite (descrip, sizeof(char), 80, pFile);
    fwrite (aux_file, sizeof(char), 24, pFile);
    fwrite (qform_code, sizeof(short), 1, pFile);
    fwrite (sform_code, sizeof(short), 1, pFile);
    fwrite (quatern_b, sizeof(float), 1, pFile);
    fwrite (quatern_c, sizeof(float), 1, pFile);
    fwrite (quatern_d, sizeof(float), 1, pFile);
    fwrite (qoffset_x, sizeof(float), 1, pFile);
    fwrite (qoffset_y, sizeof(float), 1, pFile);
    fwrite (qoffset_z, sizeof(float), 1, pFile);
    fwrite (srow_x, sizeof(float), 4, pFile);
    fwrite (srow_y, sizeof(float), 4, pFile);
    fwrite (srow_z, sizeof(float), 4, pFile);
    fwrite (intent_name, sizeof(char), 16, pFile);
    fwrite (magic, sizeof(char), 4, pFile);
    
    fseek( pFile , vox_offset[0], SEEK_SET);
    
    fwrite (odf_whole , sizeof(float) , xdim*ydim*zdim*half_no_odf_vertices , pFile);
    fclose (pFile);
    
    
    fftw_free(in_fft);
    fftw_free(out_fft);
    
    //end timer
    time_t rawtime2;
    struct tm * timeinfo2;
    time ( &rawtime2 );
    timeinfo2 = localtime ( &rawtime2 );
    cout << "END TIME/DATE: " << asctime (timeinfo2) << endl;
    
    cout << "FINISHED." << endl;
    
    return 0;
}