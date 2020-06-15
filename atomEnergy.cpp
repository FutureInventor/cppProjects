#include <iostream>
#include <vector> 
#include <stdlib.h> 
#include <random>
#include <fstream>
#include<cmath>
using namespace std;

int main() {

	int Nr = 10;
	int N = pow(Nr, 3);
	int Neq = 10000;
	double N_av = 6.022045e23;
	double m_at = 40.0;
	double m_at_kg = m_at / N_av * 1e-3;
	double ro = 1395.0;
	double L = pow((N * m_at_kg / ro), (1.0 / 3.0));
	double kb = 1.38065e-23;
	double sig = 3.41e-10;
	double eps_K = 119.8;
	double T_rz = 87.3; //temperatura zadana
	double eps_J = eps_K * kb;
	double T_zr = T_rz / eps_K; //temperatura zredukowana
	double dt = 1.0e-14;
	double N_d = N;
	double Ec;
	double sig_E;
	double T_ch_K;


	double x[1000];
	double y[1000];
	double z[1000];



	double x0 = -1.0;
	double y0 = -1.0;
	double z0 = -1.0;

	double vx[1000];
	double vy[1000];
	double vz[1000];

	for (int i = 1; i <= 1000; i++)
	{
		x[i] = 0.;
		y[i] = 0.;
		z[i] = 0.;
		vx[i] = 0.;
		vy[i] = 0.;
		vz[i] = 0.;
	}

	double jd = (1.0 / 2.0) * L;
	double sig_zr = sig / jd;
	double sig_zr_kw = sig_zr * sig_zr;

	double dt_zr = dt / sqrt(m_at_kg * jd * jd / eps_J);
	double dt_zr_kw = dt_zr * dt_zr;
	double v_max = dt_zr * sqrt(3.0 * T_zr);

	double A = 2.0 / Nr;
	double B = A / 10.0;

	int ijk = 0.0;

	for (int i = 1; i <= Nr; i++)
	{

		for (int j = 1; j <= Nr; j++)
		{

			for (int k = 1; k <= Nr; k++)
			{

				ijk++;
				random_device rd;
				mt19937 gen(rd());
				uniform_real_distribution<> dis(-0.5, 0.5);
				double r = dis(gen);
				x[ijk] = x0 + (i)*A + r * B;
				r = dis(gen);
				y[ijk] = y0 + (j)*A + r * B;
				r = dis(gen);
				z[ijk] = z0 + (k)*A + r * B;
				x[ijk] = x[ijk] - 2. * int(x[ijk]);
				y[ijk] = y[ijk] - 2. * int(y[ijk]);
				z[ijk] = z[ijk] - 2. * int(z[ijk]);
				r = dis(gen);
				vx[ijk] = 2.0 * v_max * r;
				r = dis(gen);
				vy[ijk] = 2.0 * v_max * r;
				r = dis(gen);
				vz[ijk] = 2.0 * v_max * r;
			}
		}
	}

	double vxc = 0.;
	double vyc = 0.;
	double vzc = 0.;
	for (int i = 1; i <= N; i++)
	{
		vxc = vxc + vx[i];
		vyc = vyc + vy[i];
		vzc = vzc + vz[i];
	}
	vxc = vxc / N_d;
	vyc = vyc / N_d;
	vzc = vzc / N_d;
	double Ek_zr;
	double T_ch;
	double T_ch_zr;
	Ek_zr = 0.;
	for (int i = 1; i <= N; i++)
	{
		vx[i] = vx[i] - vxc;
		vy[i] = vy[i] - vyc;
		vz[i] = vz[i] - vzc;
		Ek_zr += (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
	}
	Ek_zr = Ek_zr / dt_zr_kw;
	T_ch_zr = Ek_zr / (3.0 * (N_d - 1.0));
	T_ch_K = T_ch_zr * eps_K;

	double skal;
	if (abs(T_ch_K - T_rz) > 1.)
	{
		skal = sqrt(T_rz / T_ch_K);
		for (int i = 1; i <= N; i++)
		{
			vx[i] = vx[i] * skal;
			vy[i] = vy[i] * skal;
			vz[i] = vz[i] * skal;
		}
	}

	double E_p;
	double Fx[1000];
	double Fy[1000];
	double Fz[1000];
	double xij;
	double yij;
	double zij;
	double rij2;
	double sr2;
	double sr6;
	double dF;
	for (int i = 1; i <= N; i++)
	{
		Fx[i] = 0.0;
		Fy[i] = 0.0;
		Fz[i] = 0.0;
	}
	E_p = 0.;
	for (int i = 1; i < N; i++)
	{
		x0 = x[i];
		y0 = y[i];
		z0 = z[i];

		for (int j = i + 1; j <= N; j++)
		{
			xij = x0 - x[j];
			yij = y0 - y[j];
			zij = z0 - z[j];
			xij = xij - 2. * int(xij);
			yij = yij - 2. * int(yij);
			zij = zij - 2. * int(zij);

			rij2 = xij * xij + yij * yij + zij * zij;
			if (rij2 <= 1.0)
			{
				sr2 = sig_zr_kw / rij2;
				sr6 = sr2 * sr2 * sr2;
				E_p = E_p + sr6 * (sr6 - 1.0);
				dF = sr6 * (2.0 * sr6 - 1.0) / rij2;
				Fx[i] = Fx[i] + dF * xij;
				Fx[j] = Fx[j] - dF * xij;
				Fy[i] = Fy[i] + dF * yij;
				Fy[j] = Fy[j] - dF * yij;
				Fz[i] = Fz[i] + dF * zij;
				Fz[j] = Fz[j] - dF * zij;
			}
		}
	}

	E_p = E_p * 4.;

	for (int i = 1; i <= N; i++)
	{
		Fx[i] = Fx[i] * 24.;
		Fy[i] = Fy[i] * 24.;
		Fz[i] = Fz[i] * 24.;
	}

	double E1 = 0.;
	double E2 = 0.;

	for (int k = 1; k <= Neq; k++)
	{

		for (int i = 1; i <= N; i++)
		{
			vx[i] = vx[i] + 0.5 * Fx[i] * dt_zr_kw;
			vy[i] = vy[i] + 0.5 * Fy[i] * dt_zr_kw;
			vz[i] = vz[i] + 0.5 * Fz[i] * dt_zr_kw;
			x[i] = x[i] + vx[i];
			y[i] = y[i] + vy[i];
			z[i] = z[i] + vz[i];
			x[i] = x[i] - 2.0 * int(x[i]);
			y[i] = y[i] - 2.0 * int(y[i]);
			z[i] = z[i] - 2.0 * int(z[i]);
		}
		E_p = 0.;
		for (int i = 1; i <= N; i++)
		{
			Fx[i] = 0.0;
			Fy[i] = 0.0;
			Fz[i] = 0.0;
		}

		for (int i = 1; i < N; i++)
		{
			x0 = x[i];
			y0 = y[i];
			z0 = z[i];

			for (int j = i + 1; j <= N; j++)
			{
				xij = x0 - x[j];
				yij = y0 - y[j];
				zij = z0 - z[j];
				xij = xij - 2. * int(xij);
				yij = yij - 2. * int(yij);
				zij = zij - 2. * int(zij);

				rij2 = xij * xij + yij * yij + zij * zij;
				if (rij2 <= 1.0)
				{
					sr2 = sig_zr_kw / rij2;
					sr6 = sr2 * sr2 * sr2;
					E_p = E_p + sr6 * (sr6 - 1.0);
					dF = sr6 * (2.0 * sr6 - 1.0) / rij2;
					Fx[i] = Fx[i] + dF * xij;
					Fx[j] = Fx[j] - dF * xij;
					Fy[i] = Fy[i] + dF * yij;
					Fy[j] = Fy[j] - dF * yij;
					Fz[i] = Fz[i] + dF * zij;
					Fz[j] = Fz[j] - dF * zij;

				}
			}
		}

		E_p = E_p * 4.;

		for (int i = 1; i <= N; i++)
		{
			Fx[i] = Fx[i] * 24.;
			Fy[i] = Fy[i] * 24.;
			Fz[i] = Fz[i] * 24.;
		}

		for (int i = 1; i <= N; i++)
		{
			vx[i] = vx[i] + 0.5 * Fx[i] * dt_zr_kw;
			vy[i] = vy[i] + 0.5 * Fy[i] * dt_zr_kw;
			vz[i] = vz[i] + 0.5 * Fz[i] * dt_zr_kw;

		}

		vxc = 0.;
		vyc = 0.;
		vzc = 0.;
		for (int i = 1; i <= N; i++)
		{
			vxc = vxc + vx[i];
			vyc = vyc + vy[i];
			vzc = vzc + vz[i];
		}
		vxc = vxc / N_d;
		vyc = vyc / N_d;
		vzc = vzc / N_d;
		Ek_zr = 0.;
		for (int i = 1; i <= N; i++)
		{
			vx[i] = vx[i] - vxc;
			vy[i] = vy[i] - vyc;
			vz[i] = vz[i] - vzc;
			Ek_zr += (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
		}
		Ek_zr = Ek_zr / dt_zr_kw;
		T_ch_zr = Ek_zr / (3.0 * (N_d - 1.0));
		T_ch_K = T_ch_zr * eps_K;

		if (abs(T_ch_K - T_rz) > 1.)
		{
			skal = sqrt(T_rz / T_ch_K);
			for (int i = 1; i <= N; i++)
			{
				vx[i] = vx[i] * skal;
				vy[i] = vy[i] * skal;
				vz[i] = vz[i] * skal;
			}
		}
		Ec = 0.5 * Ek_zr + E_p;
		Ec = Ec * N_av / N_d;
		Ec = Ec * eps_J;
		double k_d = k;
		E1 = E1 + Ec;
		E2 = E2 + Ec * Ec;
		double s_E = E2 / k_d - (E1 / k_d) * (E1 / k_d);
		sig_E = sqrt(s_E);

		if (k % 100 == 0)
			cout << "k=" << k << " " << "T_ch_K=" << T_ch_K
			<< " " << "E1/k_d=" << E1 / k_d << " " << "sig_E=" << sig_E << endl;

	}

	E1 = 0.;
	E2 = 0.;

	int Ns = 5000.;

	for (int k = 1; k <= Ns; k++)
	{

		for (int i = 1; i <= N; i++)
		{
			vx[i] = vx[i] + 0.5 * Fx[i] * dt_zr_kw;
			vy[i] = vy[i] + 0.5 * Fy[i] * dt_zr_kw;
			vz[i] = vz[i] + 0.5 * Fz[i] * dt_zr_kw;
			x[i] = x[i] + vx[i];
			y[i] = y[i] + vy[i];
			z[i] = z[i] + vz[i];
			x[i] = x[i] - 2.0 * int(x[i]);
			y[i] = y[i] - 2.0 * int(y[i]);
			z[i] = z[i] - 2.0 * int(z[i]);

		}

		for (int i = 1; i <= N; i++)
		{
			Fx[i] = 0.0;
			Fy[i] = 0.0;
			Fz[i] = 0.0;
		}
		E_p = 0.;
		for (int i = 1; i < N; i++)
		{
			x0 = x[i];
			y0 = y[i];
			z0 = z[i];

			for (int j = i + 1; j <= N; j++)
			{
				xij = x0 - x[j];
				yij = y0 - y[j];
				zij = z0 - z[j];
				xij = xij - 2. * int(xij);
				yij = yij - 2. * int(yij);
				zij = zij - 2. * int(zij);

				rij2 = xij * xij + yij * yij + zij * zij;
				if (rij2 <= 1.0)
				{
					sr2 = sig_zr_kw / rij2;
					sr6 = sr2 * sr2 * sr2;
					E_p = E_p + sr6 * (sr6 - 1.0);
					dF = sr6 * (2.0 * sr6 - 1.0) / rij2;
					Fx[i] = Fx[i] + dF * xij;
					Fx[j] = Fx[j] - dF * xij;
					Fy[i] = Fy[i] + dF * yij;
					Fy[j] = Fy[j] - dF * yij;
					Fz[i] = Fz[i] + dF * zij;
					Fz[j] = Fz[j] - dF * zij;

				}
			}
		}

		E_p = E_p * 4.;

		for (int i = 1; i <= N; i++)
		{
			Fx[i] = Fx[i] * 24.;
			Fy[i] = Fy[i] * 24.;
			Fz[i] = Fz[i] * 24.;
		}

		for (int i = 1; i <= N; i++)
		{
			vx[i] = vx[i] + 0.5 * Fx[i] * dt_zr_kw;
			vy[i] = vy[i] + 0.5 * Fy[i] * dt_zr_kw;
			vz[i] = vz[i] + 0.5 * Fz[i] * dt_zr_kw;

		}

		vxc = 0.;
		vyc = 0.;
		vzc = 0.;
		for (int i = 1; i <= N; i++)
		{
			vxc = vxc + vx[i];
			vyc = vyc + vy[i];
			vzc = vzc + vz[i];
		}
		vxc = vxc / N_d;
		vyc = vyc / N_d;
		vzc = vzc / N_d;
		Ek_zr = 0.;
		for (int i = 1; i <= N; i++)
		{
			vx[i] = vx[i] - vxc;
			vy[i] = vy[i] - vyc;
			vz[i] = vz[i] - vzc;
			Ek_zr += (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
		}
		Ek_zr = Ek_zr / dt_zr_kw;
		T_ch_zr = Ek_zr / (3.0 * (N_d - 1.0));
		T_ch_K = T_ch_zr * eps_K;

		Ec = 0.5 * Ek_zr + E_p;
		Ec = Ec * N_av / N_d;
		Ec = Ec * eps_J;
		double k_d = k;
		E1 = E1 + Ec;
		E2 = E2 + Ec * Ec;
		double s_E = E2 / k_d - (E1 / k_d) * (E1 / k_d);
		sig_E = sqrt(s_E);

		if (k % 100 == 0)
			cout << "bez skalowania" << " " << "k=" << k << " " << "T_ch_K=" << T_ch_K
			<< " " << "E1/k_d=" << E1 / k_d << " " << "sig_E=" << sig_E << endl;
	}
	return 0;
}