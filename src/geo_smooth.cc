/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

/**
 * Smoothing by geometric flow.
 **/

#include <iostream>
extern "C" { 
#include  <volume_io.h>
#include <time_stamp.h> 
}

using namespace std;


// #define I(i,j,k)  (( (i) < 0 || (i) >= sizes[0] || (j) < 0 || (j) >= sizes[1] || (k) < 0 || (k) >= sizes[2] ) ? get_volume_real_min( volume ) : get_volume_real_value( volume, (i), (j), (k), 0, 0 ) )

#define I(i,j,k) ( get_volume_real_value( volume, (i), (j), (k), 0, 0 ) )

void smooth_volume( Real dt, Volume volume, Volume output )
{
    int sizes[MAX_DIMENSIONS];
    Real dx, dy, dz, dxsq, dysq, dzsq, dxdy, dxdz, dydz;
    
    {
	Real separations[MAX_DIMENSIONS];
	get_volume_separations( volume, separations );
	dx = separations[0];
	dy = separations[1];
	dz = separations[2];
    }
    dxsq = dx * dx;
    dysq = dy * dy;
    dzsq = dz * dz;
    dxdy = dx * dy;
    dxdz = dx * dz;
    dydz = dy * dz;
    
    Real Dx, Dy, Dz, Dxx, Dyy, Dzz, Dxy, Dxz, Dyz; 

    Real val_min = get_volume_real_min( volume );
    Real val_max = get_volume_real_max( volume );

    get_volume_sizes( volume, sizes );

    Real * old_val[3];
    old_val[0] = new Real[sizes[1]*sizes[2]];
    old_val[1] = new Real[sizes[1]*sizes[2]];
    old_val[2] = new Real[sizes[1]*sizes[2]];

    // pre-fetch values at i=0
    int iim1 = 0;
    int ii = 1;
    int iip1 = 2;
    int kk = 0;
    for( int j = 0;  j < sizes[1];  ++j ) {
        for( int k = 0;  k < sizes[2];  ++k ) {
            old_val[ii][kk] = get_volume_real_value( volume, 0, j, k, 0, 0 );
            kk++;
        }
    }

    for( int i = 0;  i < sizes[0];  ++i ) {

        // pre-fetch values at i+1
        kk = 0;
        if( i+1 < sizes[0] ) {
          for( int j = 0;  j < sizes[1];  ++j ) {
            for( int k = 0;  k < sizes[2];  ++k ) {
              old_val[iip1][kk] = get_volume_real_value( volume, i+1, 
                                                         j, k, 0, 0 );
              kk++;
            }
          }
        }

        for( int j = 0;  j < sizes[1];  ++j ) {
            for( int k = 0;  k < sizes[2];  ++k ) {

                int offset = j*sizes[2]+k;

		// The D's are spatial image derivatives

                if( i == 0 || i == sizes[0]-1 ) {
                  if( sizes[0] > 1 ) {
                    if( i == 0 ) {
                      Dx = ( old_val[iip1][offset] - old_val[ii][offset] ) / dx;
                    } else {
                      Dx = ( old_val[ii][offset] - old_val[iim1][offset] ) / dx;
                    }
                  } else {
                    Dx = 0.0;
                  }
                  Dxx = 0.0;
                } else {
                  Dx = ( old_val[iip1][offset] - old_val[iim1][offset] ) / 
                       ( 2.0 * dx );
                  Dxx = ( old_val[iip1][offset] -
                          2.0 * old_val[ii][offset] +
                          old_val[iim1][offset] ) / dxsq;
                }
                if( j == 0 || j == sizes[1]-1 ) {
                  if( sizes[1] > 1 ) {
                    if( j == 0 ) {
                      Dy = ( old_val[ii][(j+1)*sizes[2]+k] - 
                             old_val[ii][offset] ) / dy;
                    } else {
                      Dy = ( old_val[ii][offset] - 
                             old_val[ii][(j-1)*sizes[2]+k] ) / dy;
                    }
                  } else {
                    Dy = 0.0;
                  }
                  Dyy = 0.0;
                } else {
                  Dy = ( old_val[ii][(j+1)*sizes[2]+k] - 
                         old_val[ii][(j-1)*sizes[2]+k] ) / ( 2.0 * dy );
                  
		  Dyy = ( old_val[ii][(j+1)*sizes[2]+k] - 
                          2.0*old_val[ii][offset] +
                          old_val[ii][(j-1)*sizes[2]+k] ) / dysq;
                }
                if( k == 0 || k == sizes[2]-1 ) {
                  if( sizes[2] > 1 ) {
                    if( k == 0 ) {
                      Dz = ( old_val[ii][offset+1] - 
                             old_val[ii][offset] ) / dz;
                    } else {
                      Dz = ( old_val[ii][offset] - 
                             old_val[ii][offset-1] ) / dz;
                    }
                  } else {
                    Dz = 0.0;
                  }
                  Dzz = 0.0;
                } else {
		  Dz = ( old_val[ii][offset+1] -
                         old_val[ii][offset-1] ) / (2.0*dz);
		  Dzz = ( old_val[ii][offset+1] - 
                          2.*old_val[ii][offset] +
                          old_val[ii][offset-1] ) / dzsq;
                }
		
                if( i == 0 || i == sizes[0]-1 || j == 0 || j == sizes[1]-1 ) {
                  Dxy = 0.0;
                } else {
                  Dxy = ( old_val[iip1][(j+1)*sizes[2]+k] -
                          old_val[iip1][(j-1)*sizes[2]+k] -
                          old_val[iim1][(j+1)*sizes[2]+k] +
                          old_val[iim1][(j-1)*sizes[2]+k] ) / ( 4.0*dxdy );
                }
		
                if( i == 0 || i == sizes[0]-1 || k == 0 || k == sizes[2]-1 ) {
                  Dxz = 0.0;
                } else {
                  Dxz = ( old_val[iip1][offset+1] -
                          old_val[iip1][offset-1] -
                          old_val[iim1][offset+1] +
                          old_val[iim1][offset-1] ) / ( 4.0*dxdz );
                }
		
                if( j == 0 || j == sizes[1]-1 || k == 0 || k == sizes[2]-1 ) {
                  Dyz = 0.0;
                } else {
                  Dyz = ( old_val[ii][(j+1)*sizes[2]+k+1] -
                          old_val[ii][(j-1)*sizes[2]+k+1] -
                          old_val[ii][(j+1)*sizes[2]+k-1] +
                          old_val[ii][(j-1)*sizes[2]+k-1] ) / ( 4.0*dydz );
                }

		Real curvature_num
		    = Dxx*(Dy*Dy + Dz*Dz)
		    + Dyy*(Dx*Dx + Dz*Dz)
		    + Dzz*(Dx*Dx + Dy*Dy)
		    - 2.*(Dxy*Dx*Dy + Dxz*Dx*Dz + Dyz*Dy*Dz);

		// The denominator in the expression of curvature is the
		// gradient magnitude, cubed.  Since we ultimately multiply
		// curvature by the gradient magnitude, we compute here the
		// resulting denominator: squared gradient magnitude.
		Real den = Dx*Dx + Dy*Dy + Dz*Dz + 1.0e-10;

                Real new_val = old_val[ii][offset] + dt * curvature_num / den;
                if( new_val < val_min ) new_val = val_min;
                if( new_val > val_max ) new_val = val_max;

		set_volume_real_value( output, i, j, k, 0, 0, new_val );
	    }
	}
        iim1 = ( iim1+1 ) % 3;
        ii = ( ii+1 ) % 3;
        iip1 = ( iip1+1 ) % 3;
    }
    delete [] old_val[0];
    delete [] old_val[1];
    delete [] old_val[2];
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

    cout << N << " time steps of size " << dt << endl;
    
    Real separations[MAX_DIMENSIONS];
    get_volume_separations( in_volume, separations );
    Real dx = separations[0];
    Real dy = separations[1];
    Real dz = separations[2];

    // Check stability condition for explicit discretization.

    Real dt_max = 0.5 / ( 1.0 / ( dx * dx ) + 1.0 / ( dy * dy ) + 1.0 / ( dz * dz ) );
    if( dt >= dt_max ) {
      cout << "Time step of " << dt << " violates stability condition." << endl;
      Real tfinal = dt * N;
      N = ceil ( dt * N / dt_max );
      dt = tfinal / N;
      cout << "Adjusting to " << N << " steps of " << dt << " for stability." << endl;
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
        cerr << endl << "Copyright Alan C. Evans" << endl
                     << "Professor of Neurology" << endl
                     << "McGill University" << endl;
	return 1;
    }

    if ((ac == 6)&&(strcmp(av[5],"-outputsteps") == 0)) 
      allsteps = TRUE;

    return process_file(time_stamp(ac, av),
			atof( av[1] ), atoi( av[2] ), av[3], av[4], allsteps );
}
