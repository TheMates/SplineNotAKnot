/*
* Not-a-knot spline.
* This utilizes library Eigen, get the files from https://gitlab.com/libeigen/eigen
* It is easy to use. Add it to your project or check if the include path here coresponsd 
* to location of your files.
* 
* Inspiration from http://www.cs.tau.ac.il/~turkel/notes/numeng/spline_note.pdf
* You should check the code for MATLAB 
* Download from here http://www.cs.cornell.edu/courses/cs4210/2015fa/CVLBook/new_page_1.htm
* 
* Copyright(C) 2020 Matous Vrbik (matousvrbik[at]gmail.com)
*
*This program is free software; you can redistribute itand /or
*modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.If not, see < http://www.gnu.org/licenses/>.
*-------------------------------------------------------------------- -
*/


#ifndef SPLINE_NAK
#define SPLINE_NAK

#include <vector>
#include <cassert>


// CHECK EIGEN PATH
#include "Eigen/Eigen"

namespace SplineNaK
{
	typedef std::vector<double> Vectord;

	/**
	 * Cubic spline class, that implements not-a-knot interpolation.
	 */
	class Spline
	{
	public:
		Spline() = default;
		~Spline()= default;

		// Sets 
		void setPoints(std::vector<double>& x, std::vector<double>& y)
		{
			assert(x.size() == y.size());
			xp = &x;

			size_t n = x.size();

			Vectord Dx;						// derivative of x
			Vectord yp;						// derivative of y (dy/dx)
			Eigen::VectorXd r(n-2);			// r-side

			Eigen::SparseMatrix<double> Tsp(n - 2, n - 2);			// sparse matrix
			Tsp.reserve(Eigen::VectorXi::Constant(n - 2, 3));		//three nonzero entries in each column


			diff(x, Dx);
			diff(y, yp);
			for (int i = 0; i < yp.size(); ++i)
				yp[i] /= Dx[i];

			for(auto i = 1;i<n-3;++i)
			{

				Tsp.insert(i, i - 1) = Dx[i + 1];
				Tsp.insert(i, i) = 2 * (Dx[i] + Dx[i + 1]);
				Tsp.insert(i, i + 1) = Dx[i];

				r(i) = 3 * (Dx[i + 1] * yp[i] + Dx[i] * yp[i + 1]);
			}

			// not-a-knot computation of slopes
			double q = Dx[0] * Dx[0] / Dx[1];

			Tsp.insert(0, 0) = 2 * Dx[0] + Dx[1] + q;
			Tsp.insert(0, 1) = Dx[0] + q;


			r(0) = Dx[1] * yp[0] + Dx[0] * yp[1] + 2 * yp[1] * (q + Dx[0]);
			q = Dx[n - 2] * Dx[n-2] / Dx[n-3];
			r(n - 3) = Dx[n - 2] * yp[n - 3] + Dx[n - 3] * yp[n - 2] + 2 * yp[n - 3] * (Dx[n - 2] + q);

			Tsp.insert(n - 3, n - 3) = 2 * Dx[n - 2] + Dx[n - 3] + q;
			Tsp.insert(n - 3, n - 4) = Dx[n - 2] + q;


			Eigen::SparseLU<Eigen::SparseMatrix<double>> chol(Tsp);
			Eigen::VectorXd stilde = chol.solve(r);		//cholesky solver

			// first and last slopes
			double s0 = -stilde(0) + 2*yp[0];
			s0 = s0 + ((Dx[0]* Dx[0])/ (Dx[1] * Dx[1])) * (stilde(0) + stilde(1) - 2 * yp[1]);

			double sn = -stilde(n - 3) + 2 * yp[n - 2];
			sn = sn + ((Dx[n - 2] * Dx[n - 2]) / (Dx[n - 3] * Dx[n - 3])) * (stilde(n - 4) + stilde(n - 3) - 2 * yp[n - 3]);

			Vectord s(n);
			for (auto i = 1; i < n - 1; ++i)
				s[i] = stilde(i - 1);
			s[0] = s0;
			s[n - 1] = sn;


			// vectors of coefficients a,b,c,d
			a = Vectord(&y[0], &y[n - 1]);
			b = Vectord(&s[0], &s[n - 1]);
			c.resize(n - 1);
			d.resize(n - 1);
			for (auto i = 0; i < c.size(); ++i)
			{
				c[i] = (yp[i] - s[i]) / Dx[i];
				d[i] = (s[i + 1] + s[i] - 2 * yp[i]) / (Dx[i] * Dx[i]);
			}				
		}

		double operator ()(double x_s)
		{
			std::vector<double>::const_iterator it;
			it = std::lower_bound(xp->begin(), xp->end(), x_s);
			int idx = (std::max)(int(it - xp->begin()) - 1, 0);		//windows.h defines these macros, i want to use functions
			idx = (std::min)(idx,int(xp->size()-2));

			double y_s = d[idx] * (x_s - (*xp)[idx + 1]) + c[idx];
			y_s = y_s * (x_s - (*xp)[idx]) + b[idx];
			y_s = y_s * (x_s - (*xp)[idx]) + a[idx];

			return y_s;
		}

	private:

		// piecewise coeffs
		Vectord a;
		Vectord b;
		Vectord c;
		Vectord d;

		Vectord* xp = nullptr;	//pointer to original x data

		//support methods

		/**
		 * Calculates derivative of vector in.
		 * 
		 * out[i] = in[i+1] - in[i]
		 * 
		 * The out vector is then one shorter than in.
		 */
		static void diff(Vectord &in, Vectord& out)
		{
			out.resize(in.size() - 1);
			for (int i = 0; i < out.size(); ++i)
				out[i] = in[i + 1] - in[i];
		}
	};
}

#endif //SPLINE_NAK

