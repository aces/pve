/**
 * Smoothing by geometric flow.
 **/

#include <iostream>
extern "C" { 
#include  <volume_io.h>
#include <time_stamp.h> 
}

using namespace std;


#define I(i,j,k)  (( (i) < 0 || (i) >= sizes[0] || (j) < 0 || (j) >= sizes[1] || (k) < 0 || (k) >= sizes[2] ) ? get_volume_real_min( volume ) : get_volume_real_value( volume, (i), (j), (k), 0, 0 ) )

void smooth_volume( Real dt, Volume volume, Volume output )
{
    int sizes[MAX_DIMENSIONS];
    Real dx, dy, dz;
    
    {
	Real separations[MAX_DIMENSIONS];
	get_volume_separations( volume, separations );
	dx = separations[0];
	dy = separations[1];
	dz = separations[2];
    }
    
    get_volume_sizes( volume, sizes );
    for( int i = 0;  i < sizes[0];  ++i ) {
        for( int j = 0;  j < sizes[1];  ++j ) {
            for( int k = 0;  k < sizes[2];  ++k ) {

		// The D's are spatial image derivatives

		Real Dx = (I(i+1,j,k) - I(i-1,j,k)) / 2.*dx;
		Real Dy = (I(i,j+1,k) - I(i,j-1,k)) / 2.*dy;
		Real Dz = (I(i,j,k+1) - I(i,j,k-1)) / 2.*dz;
		
		Real Dxx = (I(i+1,j,k) - 2.*I(i,j,k) + I(i-1,j,k)) / dx*dx;
		Real Dyy = (I(i,j+1,k) - 2.*I(i,j,k) + I(i,j-1,k)) / dy*dy;
		Real Dzz = (I(i,j,k+1) - 2.*I(i,j,k) + I(i,j,k-1)) / dz*dz;

		Real Dxy = ( + I(i+1,j+1,k) - I(i+1,j-1,k)
			     - I(i-1,j+1,k) + I(i-1,j-1,k) ) / 4.*dx*dy;
		
		Real Dxz = ( + I(i+1,j,k+1) - I(i+1,j,k-1)
			     - I(i-1,j,k+1) + I(i-1,j,k-1) ) / 4.*dx*dz;
		
		Real Dyz = ( + I(i,j+1,k+1) - I(i,j-1,k+1)
			     - I(i,j+1,k-1) + I(i,j-1,k-1) ) / 4.*dy*dz;

		Real curvature_num
		    = Dxx*(Dy*Dy + Dz*Dz)
		    + Dyy*(Dx*Dx + Dz*Dz)
		    + Dzz*(Dx*Dx + Dy*Dy)
		    - 2.*(Dxy*Dx*Dy + Dxz*Dx*Dz + Dyz*Dy*Dz);

		// The denominator in the expression of curvature is the
		// gradient magnitude, cubed.  Since we ultimately multiply
		// curvature by the gradient magnitude, we compute here the
		// resulting denominator: squared gradient magnitude.
		Real den = Dx*Dx + Dy*Dy + Dz*Dz;

		set_volume_real_value( output, i, j, k, 0, 0,
				       I(i,j,k) + dt * curvature_num / den );
	    }
	}
    }
}


int process_file( char* history,
		  Real dt, int N, char* in_name, char* out_prefix, BOOLEAN allsteps ) 
{
    Volume in_volume;
    char* out_name = new char[strlen(out_prefix) + 10];

    if (allsteps) {
      cout << "Filtering from " << in_name << " to "    
	   << out_prefix << "_*.mnc" << endl;
    }

    cout << N << " time steps of size " << dt << endl;
    
    if ( input_volume( in_name, 3, NULL, 
		       MI_ORIGINAL_TYPE, 0, 0, 0,
		       TRUE, &in_volume, NULL ) != OK ) {
	return 1;
    }

    if ( get_volume_n_dimensions( in_volume ) != 3 ) {
	cerr << "Error: volume in " << in_name 
	     << " does not have three dimensions." << endl;
	return 1;
    }

    Volume out_volume 
	= copy_volume_definition( in_volume, MI_ORIGINAL_TYPE, 0, 0, 0 );

    int forward = 1;

    if (!allsteps)
    	  sprintf( out_name, "%s", out_prefix );
    

    for ( int i = 1; i <= N; ++i ) {
	cout << "Iteration " << i << endl;

	if (allsteps)
	  sprintf( out_name, "%s_%d.mnc", out_prefix, i );

        int rv = OK;
	if ( forward ) {
	    smooth_volume( dt, in_volume, out_volume );
	    if ((allsteps) || (i == N)) {
	      rv = output_modified_volume( out_name, 
					   MI_ORIGINAL_TYPE, 0, 0, 0,
					   out_volume, in_name, history, NULL );
	    }
	  } else {
	    smooth_volume( dt, out_volume, in_volume );
	    if ((allsteps) || (i == N)) {
	      rv = output_modified_volume( out_name, 
					   MI_ORIGINAL_TYPE, 0, 0, 0,
					   in_volume, in_name, history, NULL );
	    }
	  }
	forward = !forward;
	if ( rv != OK )
	  return 1;
    }
    return 0;
}



int  main( int ac, char* av[] )
{
  
  BOOLEAN allsteps = FALSE;

    if (( ac != 5 ) && ( ac != 6)) {
	cerr << "usage: " << av[0] 
	     << " dt N input output [-outputsteps]" << endl;
	return 1;
    }

    if ((ac == 6)&&(strcmp(av[5],"-outputsteps") == 0)) 
      allsteps = TRUE;

    return process_file(time_stamp(ac, av),
			atof( av[1] ), atoi( av[2] ), av[3], av[4], allsteps );
}
