#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;
const int N = 100;
const int m = 2;
const float EPS = 1e-5;
int a = 0;


void outf2(float P[m], ofstream& f3) {
	for (unsigned i = 0; i < 2; i++) f3 << left << setw(7) << P[i];
	f3 << "\n";
}

void outf(float P[N][m], unsigned a, ofstream& f3) {
	f3 << "\n\n������� ����� ����� : \n";
	for (unsigned i = 0; i < a; i++) outf2(P[i], f3);
}

void inpf(ifstream& f1, float P[N][m], ofstream& f3) {
	int c = 0;
	int koli = 0;
	float num;
	while (true)
	{
		if (f1.eof()) break;
		f1 >> num;

		if (c < N)
		{
			P[a][c] = num;
			koli++;
			f3 << P[a][c] << " ";
			c++;
		}

		if (f1.peek() == '\n' || f1.eof())
		{
			if (c > 1) {
				for (int k = 2; k < c; k++) {
					f3 << "������ �������� � ������: " << P[a][k] << ". ";
				}
			}
			if (c == 1) {
				f3 << "������ �������� 1 ����������";
				a--;

			}
			a++;

			if (f1.eof())
			{
				break;
			}
			c = 0;
			f3 << "\n";
		}
	}
}

float process(ofstream& f3, float P[N][m], int a, ofstream& f4) {
	float Ax, Ay, Bx, By, Cx, Cy, Dx, Dy, Ex, Ey, Fx, Fy; //����� ������ 
	float Lx, Ly, Hx, Hy; //��������������� �����
	int i, j, k, p, q, w, e; //���������� ��� ��������
	float AC, CB, AB, AL, LB, LC, CD, LD, DE, EH, DH, AF, FH, AH, EF; //����� ��������  
	float Kab, Kbc, Bab, Kcl, Bcl, Kcd, Kde, Kaf; //������������ ������
	int kolvo = 0, kl;
	float Perimetrs[20] = {0};
	float Pr = 0;
	float Prmax = -1;
	for (i = 0; i < a; i++) {
		Ax = P[i][0];
		Ay = P[i][1];
		f3 << "\n\n   ������ ���������� ����� ������  \n ";
		f3 << "����� ����� A(" << Ax << ";" << Ay << ")\n";
		for (j = 0; j < a; j++) {
			if ((P[j][0] == Ax) && (P[j][1] == Ay)) {
				f3 << "���� ����� B\n����� B (" << P[j][0] << ";" << P[j][1] << ") �� ��������. ���������� ������� � ������ A\n";
				continue;
			}
			else {
				Bx = P[j][0];
				By = P[j][1];
				f3 << "���� ����� B\n����� B (" << Bx << ";" << By << ") ��������. ��������� 2 ����� A(" << Ax << ";" << Ay << ")" << " B(" << Bx << "; " << By << ")\n";
			}
			for (k = 0; k < a; k++) {
				if (((P[k][0] == Ax) && (P[k][1] == Ay)) || (((P[k][0] == Bx) && (P[k][1] == By)))) {
					f3 << "���� ����� C\n����� C (" << P[k][0] << ";" << P[k][1] << ") �� ��������. ���������� ������� � ������������ ������ �����\n";
					continue;
				}
				else {
					Cx = P[k][0];
					Cy = P[k][1];
					AC = (sqrt(((Cx - Ax) * (Cx - Ax)) + ((Cy - Ay) * (Cy - Ay)))); 
					CB = (sqrt(((Cx - Bx) * (Cx - Bx)) + ((Cy - By) * (Cy - By))));  
					AB = (sqrt(((Bx - Ax) * (Bx - Ax)) + ((By - Ay) * (By - Ay)))); 
					if ((Ax != Bx) && (Ay != By)) { 
						Kab = ((By - Ay) / (Bx - Ax)); 
						Bab = Ay - (Kab * Ax); 
						Kbc = ((Cy - By) / Cx - Bx);  
						Kcl = (-((Bx - Ax) / (By - Ay))); 
						Bcl = Cy - (Kcl * Cx); 
						Lx = ((Bcl - Bab) / (Kab - Kcl)); 
						Ly = (Kab * Lx) + Bab; 
						AL = (sqrt(((Lx - Ax) * (Lx - Ax)) + ((Ly - Ay) * (Ly - Ay)))); 
						LB = (sqrt(((Lx - Bx) * (Lx - Bx)) + ((Ly - By) * (Ly - By))));  
					}
					if (Ax == Bx) {
						Lx = Ax;
						Ly = Cy;
						AL = (sqrt(((Lx - Ax) * (Lx - Ax)) + ((Ly - Ay) * (Ly - Ay)))); 
						LB = (sqrt(((Lx - Bx) * (Lx - Bx)) + ((Ly - By) * (Ly - By)))); 
					}
					if (Ay == By) {
						Lx = Cx;
						Ly = Ay;
						AL = (sqrt(((Lx - Ax) * (Lx - Ax)) + ((Ly - Ay) * (Ly - Ay)))); 
						LB = (sqrt(((Lx - Bx) * (Lx - Bx)) + ((Ly - By) * (Ly - By)))); 
					}
					f3 << "����� L � ������������: " << Lx << "   " << Ly << "\n"; 
					if (fabs(AB - (AC + CB)) < EPS) {
						f3 << "���� ����� C\n����� C (" << Cx << ";" << Cy << ") �� ��������. ��� ����� �� ������� AB, ������ ��������� ������\n";
						continue;
					}
					if ((Ax != Bx) && (Ay != By)) {
						if (fabs(Kab - Kbc) < EPS) { 
							f3 << "���� ����� C\n����� C (" << Cx << ";" << Cy << ") �� ��������. ��� ����� �� ������ AB, ������ ��������� ������\n";
							continue;
						}
					}
					if (Ax == Bx) {
						if (Cx == Ax) {
							f3 << "���� ����� C\n����� C (" << Cx << ";" << Cy << ") �� ��������. ��� ����� �� ������ AB, ������ ��������� ������\n";
							continue;
						}
					}
					if (Ay == By) {
						if (Cy == Ay) {
							f3 << "���� ����� C\n����� C (" << Cx << ";" << Cy << ") �� ��������. ��� ����� �� ������ AB, ������ ��������� ������\n";
							continue;
						}
					}
					if (((Lx == Ax) && (Ly == Ay)) || ((Lx == Bx) && (Ly == By)) || ((Lx == Cx) && (Ly == Cy))) {
						f3 << "���� ����� C\n����� C (" << Cx << ";" << Cy << ") �� ��������. ������ CL �������� ����� ����� A ��� ����� ����� B, ��� ����� C ��������� � ������ L ������ ��������� ������\n";
						continue;
					}
					if (fabs(AB - (AL + LB)) < EPS) {
						f3 << "���� ����� C\n����� C (" << Cx << ";" << Cy << ") �� ��������. ������ CL ���������� ������� AB, ������ ��������� ������\n";
						continue;
					}
					if (fabs(AL - (AB + LB)) > EPS) {
						f3 << "���� ����� C\n����� C (" << Cx << ";" << Cy << ") �� ��������. ����� C ����� �� � ��� �������, ������ ��������� ������\n";
						continue;
					}
					if ((((Cx - Ax) * (By - Ay)) - ((Cy - Ay) * (Bx - Ax))) < 0) {
						f3 << "���� ����� C\n����� C (" << Cx << ";" << Cy << ") �� ��������. ����� C ����� �� � ��� �������, ������ ��������� ������\n";
						continue;
					}
					else {
						f3 << "����� C(" << Cx << ";" << Cy << ") ��������. ��������� 3 ����� A(" << Ax << ";" << Ay << ")" << " B(" << Bx << "; " << By << ")" << " C(" << Cx << ";" << Cy << ")" << "\n";
						for (q = 0; q < a; q++) {
							if (((P[q][0] == Ax) && (P[q][1] == Ay)) || ((P[q][0] == Bx) && (P[q][1] == By)) || ((P[q][0] == Cx) && (P[q][1] == Cy)) || ((P[q][0] == Lx) && (P[q][1] == Ly))) {
								f3 << "���� ����� D\n����� D (" << P[q][0] << ";" << P[q][1] << ") �� ��������. ���������� ������� � ������������ ������ �����\n";
								continue;
							}
							else {
								Dx = P[q][0];
								Dy = P[q][1];
								LC = (sqrt(((Cx - Lx) * (Cx - Lx)) + ((Cy - Ly) * (Cy - Ly)))); 
								CD = (sqrt(((Cx - Dx) * (Cx - Dx)) + ((Cy - Dy) * (Cy - Dy)))); 
								LD = (sqrt(((Lx - Dx) * (Lx - Dx)) + ((Ly - Dy) * (Ly - Dy)))); 
								if (Cx != Dx) {
									Kcd = ((Dy - Cy) / (Dx - Cx));
								}
								if ((Ax != Bx) && (Ay != By)) {
									if ((fabs(Kcl - Kcd)) > EPS) { 
										f3 << "���� ����� D\n����� D (" << Dx << ";" << Dy << ") �� ��������. ����� D �� ����� �� ������ CL, ������ ��������� ������\n";
										continue;
									}
								}
								if (Ax == Bx) {
									if (Cy != Dy) {
										f3 << "���� ����� D\n����� D (" << Dx << ";" << Dy << ") �� ��������. ����� D �� ����� �� ������ CL, ������ ��������� ������\n"; 
										continue;
									}
								}
								if (Ay == By) {
									if (Cx != Dx) {
										f3 << "���� ����� D\n����� D (" << Dx << ";" << Dy << ") �� ��������. ����� D �� ����� �� ������ CL, ������ ��������� ������\n";
										continue;
									}
								}
								if ((fabs(LD - (LC + CD)) > EPS)) { 
									f3 << "���� ����� D\n����� D (" << Dx << ";" << Dy << ") �� ��������. ����� D ����� �� ���, ��� �����, ������ ��������� ������\n";
									continue;
								}
								else {
									f3 << "����� D(" << Dx << ";" << Dy << ") ��������. ��������� 4 ����� A(" << Ax << ";" << Ay << ")" << " B(" << Bx << "; " << By << ")" << " C(" << Cx << ";" << Cy << ")" << " D(" << Dx << ";" << Dy << ")" << "\n";
									for (p = 0; p < a; p++) {
										if (((P[p][0] == Ax) && (P[p][1] == Ay)) || ((P[p][0] == Bx) && (P[p][1] == By)) || ((P[p][0] == Cx) && (P[p][1] == Cy)) || ((P[p][0] == Dx) && (P[p][1] == Dy)) || ((P[p][0] == Dx) && (P[p][1] == Dy))) {
											f3 << "���� ����� E\n����� E (" << P[p][0] << ";" << P[p][1] << ") �� ��������. ���������� ������� � ������������ ������ �����\n";
											continue;
										}
										else {
											Ex = P[p][0];
											Ey = P[p][1];
											Hx = Ax + Dx - Lx;
											Hy = Ay + Dy - Ly;
											DE = (sqrt(((Dx - Ex) * (Dx - Ex)) + ((Dy - Ey) * (Dy - Ey))));
											EH = (sqrt(((Ex - Hx) * (Ex - Hx)) + ((Ey - Hy) * (Ey - Hy))));
											DH = (sqrt(((Hx - Dx) * (Hx - Dx)) + ((Hy - Dy) * (Hy - Dy))));
											if (Ex != Dx) {
												Kde = ((Ey - Dy) / (Ex - Dx));
											}
											if ((Ax != Bx) && (Ay != By)) {
												if ((fabs(Kde - Kab)) > EPS) {
													f3 << "���� ����� E\n����� E (" << Ex << ";" << Ey << ") �� ��������. ������ DE �� ����������� ������ AB, ������ ��������� ������\n";
													continue;
												}
											}
											if (Ax == Bx) {
												if (Ex != Dx) {
													f3 << "���� ����� E\n����� E (" << Ex << ";" << Ey << ") �� ��������. ������ DE �� ����������� ������ AB, ������ ��������� ������\n";
													continue;
												}
											}
											if (Ay == By) {
												if (Ey != Dy) {
													f3 << "���� ����� E\n����� E (" << Ex << ";" << Ey << ") �� ��������. ������ DE �� ����������� ������ AB, ������ ��������� ������\n";
													continue;
												}
											}
											if ((fabs(DE - AB)) > EPS) {
												f3 << "���� ����� E\n����� E (" << Ex << ";" << Ey << ") �� ��������. ����� ������� DE �� ����� ����� ������� AB, ������ ��������� ������\n";
												continue;
											}
											if ((fabs(DH - (DE + EH)) > EPS)) {
												f3 << "���� ����� E\n����� E (" << Ex << ";" << Ey << ") �� ��������. ����� E ����� �� ���, ��� �����, ������ ��������� ������\n";
												continue;
											}
											else {
												f3 << "����� E(" << Ex << ";" << Ey << ") ��������. ��������� 5 ����� A(" << Ax << ";" << Ay << ")" << " B(" << Bx << "; " << By << ")" << " C(" << Cx << ";" << Cy << ")" << " D(" << Dx << ";" << Dy << ")" << " E(" << Ex << ";" << Ey << ")" << "\n";
												for (w = 0; w < a; w++) {
													if (((P[w][0] == Ax) && (P[w][1] == Ay)) || ((P[w][0] == Bx) && (P[w][1] == By)) || ((P[w][0] == Cx) && (P[w][1] == Cy)) || ((P[w][0] == Dx) && (P[w][1] == Dy)) || ((P[w][0] == Lx) && (P[w][1] == Ly)) || ((P[w][0] == Hx) && (P[w][1] == Hy)) || ((P[w][0] == Ex) && (P[w][1] == Ey))) {
														f3 << "���� ����� F\n����� F (" << P[w][0] << ";" << P[w][1] << ") �� ��������. ���������� ������� � ������������ ������ �����\n";
														continue;
													}
													else {
														Fx = P[w][0];
														Fy = P[w][1];
														AF = (sqrt(((Ax - Fx) * (Ax - Fx)) + ((Ay - Fy) * (Ay - Fy)))); 
														FH = (sqrt(((Fx - Hx) * (Fx - Hx)) + ((Fy - Hy) * (Fy - Hy)))); 
														AH = (sqrt(((Hx - Ax) * (Hx - Ax)) + ((Hy - Ay) * (Hy - Ay))));
														if (Ax != Fx) {
															Kaf = ((Ay - Fy) / (Ax - Fx)); 
														}
														if ((Ax != Bx) && (Ay != By)) {
															if ((fabs(Kcd - Kaf)) > EPS) {
																f3 << "���� ����� F\n����� F (" << Fx << ";" << Fy << ") �� ��������. ������ AF �� ����������� ������ CD, ������ ��������� ������\n";
																continue;
															}
														}
														if (Ax == Bx) {
															if (Fy != Ay) {
																f3 << "���� ����� F\n����� F (" << Fx << ";" << Fy << ") �� ��������. ������ AF �� ����������� ������ CD, ������ ��������� ������\n";
																continue;
															}
														}
														if (Ay == By) {
															if (Fx != Ax) {
																f3 << "���� ����� F\n����� F (" << Fx << ";" << Fy << ") �� ��������. ������ AF �� ����������� ������ CD, ������ ��������� ������\n";
																continue;
															}
														}
														if ((fabs(AF - CD)) > EPS) {
															f3 << "���� ����� F\n����� F (" << Fx << ";" << Fy << ") �� ��������. ����� ������� AF �� ����� ����� ������� CD, ������ ��������� ������\n";
															continue;
														}
														if (fabs(AH - (AF + FH)) > EPS) {
															f3 << "���� ����� F\n����� F (" << Fx << ";" << Fy << ") �� ��������. ����� F ����� �� ���, ��� �����, ������ ��������� ������\n";
															continue;
														}
														else {
															f3 << "����� F(" << Fx << ";" << Fy << ") ��������. ��������� 6 ����� A(" << Ax << ";" << Ay << ")" << " B(" << Bx << "; " << By << ")" << " C(" << Cx << ";" << Cy << ")" << " D(" << Dx << ";" << Dy << ")" << " E(" << Ex << ";" << Ey << ")" << " F(" << Fx << ";" << Fy << ")" << "\n"; 
															//f4 << "������� ������: " << "A(" << Ax << "; " << Ay << ")" << " B(" << Bx << "; " << By << ")" << " C(" << Cx << "; " << Cy << ")" << " D(" << Dx << "; " << Dy << ")" << " E(" << Ex << "; " << Ey << ")" << " F(" << Fx << "; " << Fy << ")" << "\n"; 
															kolvo++;
															f4<< "������������ ������ �����: " << kolvo << "   " << "A(" << Ax << "; " << Ay << ")" << " B(" << Bx << "; " << By << ")" << " C(" << Cx << "; " << Cy << ")" << " D(" << Dx << "; " << Dy << ")" << " E(" << Ex << "; " << Ey << ")" << " F(" << Fx << "; " << Fy << ")" << "\n";
															EF = (sqrt(((Ex - Fx) * (Ex - Fx)) + ((Ey - Fy) * (Ey - Fy))));
															Pr = AB + CB + CD + DE + EF + AF;
															if (Pr > Prmax) {
																Prmax = Pr;
																kl = kolvo;
															}
															Perimetrs[kolvo] = Pr;
															f3 << "\n�������� ������ �����: " << Pr <<"\n\n";
															f4 << "�������� ������ ����� " << kolvo << " �����: " << Pr << "\n\n";
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	if (Prmax == -1) {
		f4 << "������ �� �������" << "\n";
	}
	else {
		f4 << "������������ �������� Pr = " << Prmax << " ����� ������ � ��������: ";
		for (e = 0; e <= kolvo + 1; e++) {
			if (Perimetrs[e] == Prmax) {
				f4 << e << "  ";
			}
		}
	}
	return 0;
}

int main()
{
	setlocale(LC_ALL, "Russian");
	ifstream f1("in.txt", ios::in);
	if (!f1.is_open())
	{
		cout << "������ �������� ��������� �����" << std::endl;
		return 0;
	}
	float P[N][m];

	ofstream f3("protocol.txt", ios::out);
	if (!f3.is_open())
	{
		cout << "������ �������� ����� ������" << "\n";
		return 0;
	}
	ofstream f4("out.txt", ios::out);
	if (!f4.is_open())
	{
		cout << "������ �������� ����� ������" << "\n";
		return 0;
	}

	inpf(f1, P, f3);
	f1.close();
	outf(P, a, f3);
	process(f3, P, a, f4);
	f3.close();
	return 0;
}