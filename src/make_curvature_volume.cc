/**
 * Compute curvature from image intensities.
 * Uses method of Thirion and Gourdon (INRIA Tech Report 1881-1).
 **/

#include <iostream>
#include <math.h>


extern "C" { 
#   include  <volume_io.h> 
#include <time_stamp.h>
}

const int continuity = 2;

#define PRINT(x) std::cout << #x << " = " << x << std::endl



/*! \brief Compute principal curvatures at each voxel.
 *
 * \param vol_in the input volume
 * \param vol_k1 output volume containing the maximal curvature
 * \param vol_k2 output volume containing the minimal curvature
 *
 * For \a vol_k1 and \a vol_k2, maximal and minimal refer to the
 * absolute value.  So abs(k1) >= abs(k2) at each voxel.
 *
 * To get the Gaussian curvature, simply multiply \a vol_k1
 * and \a vol_k2 voxel-wise.  To get the mean curvature,
 * form the average (k1+k2)/2.
 */
void compute_curvatures( Volume vol_in, 
			 Volume vol_k1, Volume vol_k2 )
{
    int sizes[MAX_DIMENSIONS];
    Real voxel[MAX_DIMENSIONS];
    Real value;
    Real** first_deriv;
    Real*** second_deriv;

    get_volume_sizes( vol_in, sizes );

    ALLOC2D( first_deriv, 1, N_DIMENSIONS );
    ALLOC3D( second_deriv, 1, N_DIMENSIONS, N_DIMENSIONS );

    // Track the data range for the two output volumes.
    Real min_k1, max_k1, min_k2, max_k2;

    // Visual progress bar
    progress_struct  progress;
    int steps_completed = 0;
    initialize_progress_report( &progress, FALSE,
                                sizes[0]*sizes[1]*sizes[2],
                                "Creating Curvatures" );


    for( int i = 0;  i < sizes[0];  ++i ) 
    {
	voxel[0] = i;
        for( int j = 0;  j < sizes[1];  ++j ) 
	{
	    voxel[1] = j;
            for( int k = 0;  k < sizes[2];  ++k ) 
	    {
		voxel[2] = k;
		
		evaluate_volume( vol_in, voxel, NULL, continuity, FALSE,
				 0, &value, first_deriv, second_deriv );

		// TODO: verify the following mappings.

		Real fx = first_deriv[0][0];
		Real fy = first_deriv[0][1];
		Real fz = first_deriv[0][2];

		Real fxx = second_deriv[0][0][0];
		Real fxy = second_deriv[0][0][1];
		Real fxz = second_deriv[0][0][2];

		//Real fyx = second_deriv[0][1][0];
		Real fyy = second_deriv[0][1][1];
		Real fyz = second_deriv[0][1][2];

		//Real fzx = second_deriv[0][2][0];
		//Real fzy = second_deriv[0][2][1];
		Real fzz = second_deriv[0][2][2];

		// E = e/fz^2;
		// F = f/fz^2;
		// G = g/fz^2;
		
		Real e = fx*fx + fz*fz;
		Real f = fx*fy;
		Real g = fy*fy + fz*fz;

		// L = l / sqrt(H)*fz^3;
		// M = m / sqrt(H)*fz^3;
		// N = n / sqrt(H)*fz^3;

		Real l = 2.*fx*fz*fxz - fx*fx*fzz - fz*fz*fxx;
		Real m = fx*fz*fyz + fy*fz*fxz - fx*fy*fzz - fz*fz*fxy;
		Real n = 2.*fy*fz*fyz - fy*fy*fzz - fz*fz*fyy;

		/*
		// The matrix /a b\         /E F\-1 /L M\ 
		//            |   | = scale |   |   |   | 
		//            \c d/         \F G/   \M N/
		*/

		Real scale = pow(fx*fx + fy*fy + fz*fz, 1.5) * fz*fz;

		Real a = g*l - f*m;
		Real b = g*m - f*n;
		Real c = e*m - f*l;
		Real d = e*n - f*m;

		/*
		// The principal curvatures are given by the eigenvalues of
		//       1   /a b\
		// W = ----- |   |
		//     scale \c d/
		//
		// The eigenvalues are expressible in terms of 
		// half the trace (mean curvature) and the determinant
		// (Gaussian curvature) of W.
		*/

		// These are unscaled ...
		Real Kmean = 0.5 * (a+d);
		Real Kgauss = a*d - b*c;
		
		Real delta = sqrt(Kmean*Kmean - Kgauss);
		
		Real k1 = (Kmean + delta)/scale;
		Real k2 = (Kmean - delta)/scale;

		// Gross hack!
		/* The formulae for curvature were derived by assuming
		 * a parameterization of the (iso)surface of the form z(x,y).  
		 * This assumption doesn't hold when the isosurface is vertical
		 * (hence fz = 0).  In the code, this results in a division by
		 * scale = 0.  We hide the problem by setting curvatures to 
		 * zero when the scale is small.  This is not a cure, but
		 * it does allow us to visualize the values, because there
		 * won't be infinities in the volume.
		 */
		if ( fabs(scale) < 1e-8 ) {
#if 0
		    PRINT(fx);
		    PRINT(fy);
		    PRINT(fz);
		    PRINT(scale);

		    PRINT(delta);
		    PRINT(Kmean+delta);
		    PRINT(Kmean-delta);

		    PRINT(k1);
		    PRINT(k2);
#endif
		    k1 = k2 = 0;
		}

		if ( Kmean < 0 ) {
		    Real tmp = k1;
		    k1 = k2;
		    k2 = tmp;
		}

		set_volume_real_value( vol_k1, i,j,k,0,0, k1 );
		set_volume_real_value( vol_k2, i,j,k,0,0, k2 );

		if ( i == 0 && j == 0 && k == 0 ) {
		    min_k1 = max_k1 = k1;
		    min_k2 = max_k2 = k2;
		} else {
		    if ( k1 < min_k1 ) 
			min_k1 = k1;
		    else if ( k1 > max_k1 )
			max_k1 = k1;
		    if ( k2 < min_k2 )
			min_k2 = k2;
		    else if ( k2 > max_k2 )
			max_k2 = k2;
		}

		update_progress_report( &progress, ++steps_completed );
	    }
	}
    }
    terminate_progress_report( &progress );

    set_volume_real_range( vol_k1, min_k1-1.0, max_k1+1.0 );
    set_volume_real_range( vol_k2, min_k2-1.0, max_k2+1.0 );

    FREE2D( first_deriv );
    FREE3D( second_deriv );
}


int  main( int ac, char* av[] )
{
    if ( ac != 3 ) {
	std::cerr << "usage: " << av[0] 
		  << " input out_prefix" << std::endl;
	return 1;
    }

    char* history = time_stamp( ac, av);

    STRING name_in = av[1];
    STRING out_prefix = av[2];

    STRING name_k1 = concat_strings( out_prefix, "_k1.mnc" );
    STRING name_k2 = concat_strings( out_prefix, "_k2.mnc" );

    Volume vol_in;
    
    if ( input_volume( name_in, 3, NULL, 
		       MI_ORIGINAL_TYPE, 0, 0, 0,
		       TRUE, &vol_in, NULL ) != OK ) {
	std::cerr << "error reading input volume: " << name_in << std::endl;
	return 1;
    }

    if ( get_volume_n_dimensions( vol_in ) != 3 ) {
	std::cerr << "Error: volume in " << name_in 
		  << " does not have three dimensions." << std::endl;
	return 1;
    }

    Volume vol_k1 = copy_volume_definition( vol_in, NC_FLOAT, 0, 0, 0 );
    Volume vol_k2 = copy_volume_definition( vol_in, NC_FLOAT, 0, 0, 0 );

    compute_curvatures( vol_in, vol_k1, vol_k2 );


    if ( output_modified_volume( name_k1, 
				 MI_ORIGINAL_TYPE, 0, 0, 0,
				 vol_k1, name_in, history, NULL ) != OK ) {
	std::cerr << "Output failed for volume " << name_k1 << std::endl;
    }

    if ( output_modified_volume( name_k2, 
				 MI_ORIGINAL_TYPE, 0, 0, 0,
				 vol_k2, name_in, history, NULL ) != OK ) {
	std::cerr << "Output failed for volume " << name_k2 << std::endl;
    }

    return 0;
}
