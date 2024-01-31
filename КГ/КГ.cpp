#include <stdlib.h>
#include "glut.h"
#include <vector>
#include <iostream>
#include <locale.h>
#include <string>
#include <stdio.h>
#include <string.h>

using namespace std;

GLint Width = 1200, Height = 900;
GLubyte PointSize = 5;
double m; // размерность глобальной матрицы

double n = Width * 10; // размер сетки
double h = 20;// шаг по сетке
double scale = 1;// коэффициент масштаба 

int k = 1; // количество конечных элементов
vector<vector<double>> M; // матрица для СЛАУ
vector<double> b; // вектор правой части
vector<double> q; // искомое решение
vector<vector<double>> nvtr; //nvtr[0] - левая граница кэ; nvtr[1] - правая граница кэ    

double alpha = 1;
double beta = 1;
double poz_x = 0;
double poz_y = 0;

enum keys {
	Empty,
	KeyR, KeyG, KeyB,
	KeyW, KeyA, KeyS, KeyD,
	KeyU, KeyI,	KeyGLB,
	KeyS1, KeyS2, KeyS3, KeyS4,
	Key1, Key2, Key3, Key4,
	KeyPr
};


GLubyte ColorR = 0, ColorG = 255, ColorB = 255; // цвет точки

GLubyte ColorRS = 247, ColorGS = 97, ColorBS = 79; // цвет сплайна


const double Malpha[4][4] = 
{
	{36, 3,  -36, 3}, 
	{3, 4, -3, -1},
	{-36, -3, 36, -3},	 
	{3, -1, -3, 16}
};

const double Mbeta[4][4] = 
{
	{60, 30, -60, 30},  
	{30, 16, -30, 14},
	{-60, -30, 60,-30}, 
	{30, 14, -30, 16}
};


// структура точка
struct Point
{
	double x, y; // координаты точки
};
Point tmpPoint;
vector <Point> Points; // траектория
vector<Point> pointsSpline; // сплайн

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//базисные функции
double phi(double x, int i)
{
	if (i == 0)
		return 1 - (3 * (x * x)) + (2 * (x * x * x));
	if (i == 1)
		return x - (2 * (x * x)) + (x * x * x);
	if (i == 2)
		return (3 * (x * x)) - (2 * (x * x * x));
	if (i == 3)
		return -(x * x) + (x * x * x);
}
// xx - координата, k - номер КЭ, i - номер базисной функции
double psi(double xx, int k, int i)
{
	int k1 = nvtr[0][k];
	int k2 = nvtr[1][k];
	double h = abs(Points[k2].x - Points[k1].x);
	double ksi = (xx - Points[k1].x) / h;
	return phi(ksi, i);
}
//вклад от слагаемого с alpha 
double Malphak(int i, int j, int kk)
{
	int k1 = nvtr[0][kk];
	int k2 = nvtr[1][kk];
	double h = abs(Points[k2].x - Points[k1].x);
	double q = ((i == 0 || i == 2) && (j == 0 || j == 2)) ? 1 / h :
		((i == 1 || i == 3) && (j == 1 || j == 3)) ? h :
		1;
	return Malpha[i][j] * q / 30;
}
//вклад от слагаемого с beta(формула 4.110)
double  Mbetak(int i, int j, int kk)
{
	int k1 = nvtr[0][kk];
	int k2 = nvtr[1][kk];
	double h = abs(Points[k2].x - Points[k1].x);
	double q = ((i == 0 || i == 2) && (j == 0 || j == 2)) ? 1 / (h * h * h) :
		((i == 1 || i == 3) && (j == 1 || j == 3)) ? 1 / h :
		1 / (h * h);
	return Mbeta[i][j] * q;

}

// kk - номер текущего КЭ, где находится эта точка
// o - номер этой точки
// расчет локальной матрицы
double m_loc(int i, int j, int kk)
{
	double Aij = 0;
	double xo = Points[nvtr[0][kk]].x; // последний узел не учитываем
	for (int s = nvtr[0][kk],
		oe = ((k - 1) != kk) ? (nvtr[1][kk] - 1) : nvtr[1][kk];
		s <= oe; s++, xo++)
		Aij += psi(xo, kk, i) * psi(xo, kk, j);
	return Aij + alpha * Malphak(i, j, kk) + beta * Mbetak(i, j, kk);
}
// локальный вектор правой части
double b_loc(int i, int kk)
{
	double  Fi = 0;
	double xo = Points[nvtr[0][kk]].x; // последний узел не учитываем

	for (int j = nvtr[0][kk], oe = ((k - 1) != kk) ? (nvtr[1][kk] - 1) : nvtr[1][kk];
		j <= oe; j++, xo++)
		Fi += psi(xo, kk, i) * Points[j].y;

	return Fi;
}
// глобальная матрица
void globalMatrixComplite()
{
	for (int i = 0; i < m; i++)
		for (int j = 0; j < m; j++)
			M[i][j] = 0;

	for (int kk = 0; kk < k; kk++) 
	{
		int ij0 = kk * 2;
		for (int ii = 0, i = ij0; ii < 4; ii++, i++)
			for (int jj = 0, j = ij0; jj < 4; jj++, j++)
				M[i][j] += m_loc(ii, jj, kk);
	}
}
// глобальный вектор правой части
void globalVectorComplite()
{
	for (int i = 0; i < m; i++)
		b[i] = 0;

	for (int i = 0; i < k; i++) 
	{
		int i0 = i * 2;
		for (int j = 0; j < 4; j++, i0++)
			b[i0] += b_loc(j, i);
	}
}

// получение вектора
vector<Point> getVectorPointsOutSpline()
{
	vector<Point> pointsSpline;

	double f;
	for (int i = 0; i < k; i++) 
	{
		int i1 = nvtr[0][i];
		int i2 = nvtr[1][i];
		for (int j = i1; j < i2; j++) 
		{
			double x1 = Points[j].x;
			double x2 = Points[j + 1].x;
			double h = (x2 - x1) / 20;
			for (int j = 0; j < 20; j++, x1 += h)
			{
				int p = i * 2;
				f = 0;
				for (int ii = 0; ii < 4; ii++, p++)
					f += q[p] * psi(x1, i, ii);
				Point a;
				a.x = x1;
				a.y = f;
				pointsSpline.push_back(a);
			}
		}
	}
	return pointsSpline;
}
void spline()
{
	nvtr.resize(2);
	nvtr[0].resize(k);
	nvtr[1].resize(k);
	nvtr[0][0] = 0;
	double interval = Points.size() / k;
	for (int i = 0; i < k - 1; i++) 
	{
		nvtr[1][i] = interval * (i + 1);
		nvtr[0][i + 1] = nvtr[1][i];
	}
	nvtr[1][k - 1] = Points.size() - 1;

	M.resize(m);
	b.resize(m);
	q.resize(m);
	for(int i = 0; i < m; i++)
		M[i].resize(m);

	globalMatrixComplite();
	globalVectorComplite();
}

// Гаусс
void permuteLine(int x, int y)
{
	for (int i = 0; i < m; i++)
	{
		double bufM = M[x][i];
		M[x][i] = M[y][i];
		M[y][i] = bufM;
	}
	double bufF = b[x];
	b[x] = b[y];
	b[y] = bufF;
}

int searchMinInFirstElemColumn(int c, int ls)
{
	double min = M[ls][c];
	int ind = ls;
	for (int i = ls; i < m - 1; i++)
		if (abs(min) > abs(M[i + 1][c]) && min != 0 && M[i + 1][c] != 0)
		{
			min = M[i + 1][c];
			ind = i + 1;
		}
	if (min == 0)
		return -1;
	else
		return ind;
	return ls;
}
void adductionTriangleView()
{
	for (int i = 0; i < m - 1; i++)
	{
		int MinLine = searchMinInFirstElemColumn(i, i);
		if (MinLine != -1)
		{
			if (MinLine != i)
				permuteLine(i, MinLine);
			for (int j = i + 1; j < m; j++)
				if (M[j][i] != 0) 
				{
					double Factor = -M[j][i] / M[i][i];
					b[j] += Factor * b[i];
					for (int k = i; k < m; k++)
						M[j][k] += Factor * M[i][k];
				}
		}
		else break;
	}
}
void backStroke()
{
	for (int j = m - 1; j >= 0; j--)
	{
		q[j] = b[j] / M[j][j];
		for (int i = j - 1; i >= 0; i--)
			b[i] -= M[i][j] * q[j];
	}
}
void gauss()
{
	adductionTriangleView();
	backStroke();
}

// рисование сетки
void drawGrid()
{
	glLineWidth(1);
	glColor3ub(200, 200, 200);
	glBegin(GL_LINES);
	for (double i = 0; i < n; i += h)
	{
		// сетка OX
		glVertex2f(-Width * 2 + poz_x, Height / 2 + poz_y + i);
		glVertex2f(2 * Width + poz_x, Height / 2 + i + poz_y);

		glVertex2f(-Width * 2 + poz_x, Height / 2 - i + poz_y);
		glVertex2f(2 * Width + poz_x, Height / 2 - i + poz_y);
		// сетка OY
		glVertex2f(Width / 2 + i + poz_x, -Height * 2 + poz_y);
		glVertex2f(Width / 2 + i + poz_x, 2 * Height + poz_y);

		glVertex2f(Width / 2 - i + poz_x, -Height * 2 + poz_y);
		glVertex2f(Width / 2 - i + poz_x, 2 * Height + poz_y);
	}
	glEnd();

	glLineWidth(3);
	glColor3ub(0, 0, 255);
	glBegin(GL_LINES);
	// ось OX
	glVertex2f(-Width * 2 + poz_x, Height / 2 + poz_y);
	glVertex2f(2 * Width + poz_x, Height / 2 + poz_y);
	// ось OY
	glVertex2f(Width / 2 + poz_x, -Height * 2 + poz_y);
	glVertex2f(Width / 2 + poz_x, 2 * Height + poz_y);
	glEnd();
	// стрелки
	glBegin(GL_TRIANGLES);
	// OX
	glVertex2f(Width, Height / 2 + poz_y);
	glVertex2f(Width - 20, Height / 2 + 7 + poz_y);
	glVertex2f(Width - 20 , Height / 2 - 7 + poz_y);
	// OY    
	glVertex2f(Width / 2 + poz_x, Height);
	glVertex2f(Width / 2 - 7 + poz_x, Height - 20);
	glVertex2f(Width / 2 + 7 + poz_x, Height - 20);
	glEnd();
}

// номера координат
void drawGridNumber()
{
	char buffer[_CVTBUFSIZE];
	// ось OX
	for (int i = 0; i <= Width; i += 2 * h)
	{
		// > 0
		glColor3ub(0, 0, ColorB);
		glRasterPos2f(Width / 2 + i + poz_x, Height / 2 - 15 + poz_y);
		_gcvt_s(buffer, _CVTBUFSIZE, i * scale, 5);
		for (unsigned int j = 0; buffer[j] != ',' && buffer[0] != '0'; ++j)
		{
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, buffer[j]);
		}
		glLineWidth(3);
		glColor3ub(0, 0, 255);
		glBegin(GL_LINES);
		glVertex2f(Width / 2 + i + poz_x, Height / 2 - 3 + poz_y);
		glVertex2f(Width / 2 + i + poz_x, Height / 2 + 3 + poz_y);
		glEnd();
		// < 0
		glColor3ub(0, 0, ColorB);
		glRasterPos2f(Width / 2 - i + poz_x, Height / 2 - 15 + poz_y);
		_gcvt_s(buffer, _CVTBUFSIZE, -i * scale, 5);
		for (unsigned int j = 0; buffer[j] != ',' && buffer[0] != '0'; ++j)
		{
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, buffer[j]);
		}
		glLineWidth(3);
		glColor3ub(0, 0, 255);
		glBegin(GL_LINES);
		glVertex2f(Width / 2 - i + poz_x, Height / 2 - 3 + poz_y);
		glVertex2f(Width / 2 - i + poz_x, Height / 2 + 3 + poz_y);
		glEnd();
	}

	// ось OY
	for (int i = 0; i <= Height; i += 2 * h)
	{
		// >0
		glColor3ub(0, 0, ColorB);
		glRasterPos2f(Width / 2 + 5 + poz_x, Height / 2 + i + poz_y);
		_gcvt_s(buffer, _CVTBUFSIZE, i * scale, 5);
		for (unsigned int j = 0; buffer[j] != ','; ++j)
		{
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, buffer[j]);
		}
		glLineWidth(3);
		glColor3ub(0, 0, 255);
		glBegin(GL_LINES);
		glVertex2f(Width / 2 + 3 + poz_x, Height / 2 + i + poz_y);
		glVertex2f(Width / 2 - 3 + poz_x, Height / 2 + i + poz_y);
		glEnd();
		// <0
		glColor3ub(0, 0, ColorB);
		glRasterPos2f(Width / 2 + 5 + poz_x, Height / 2 - i + poz_y);
		_gcvt_s(buffer, _CVTBUFSIZE, -i * scale, 5);
		for (unsigned int j = 0; buffer[j] != ','; ++j)
		{
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, buffer[j]);
		}
		glLineWidth(3);
		glColor3ub(0, 0, 255);
		glBegin(GL_LINES);
		glVertex2f(Width / 2 + 3 + poz_x, Height / 2 - i + poz_y);
		glVertex2f(Width / 2 - 3 + poz_x, Height / 2 - i + poz_y);
		glEnd();
	}
}

// рисование точек и траектории
void drawPoint()
{
	glColor3ub(ColorR, ColorG, ColorB);
	glPointSize(5);

	glBegin(GL_POINTS);
	for (int i = 0; i < Points.size(); i++)
		glVertex2i(Points[i].x, Points[i].y);
	glEnd();

	glLineWidth(1);
	glEnable(GL_LINE_STIPPLE);
	glLineStipple(1, 0x00FF);
	glBegin(GL_LINE_STRIP); //рисуем замкнутую линию 
	for (int i = 0; i < Points.size() && Points.size() > 1; i++)
		glVertex2i(Points[i].x, Points[i].y);
	glEnd();
	glDisable(GL_LINE_STIPPLE);
}

// рисование сплайна
void drawSpline()
{
	if (Points.size() >= 4)
	{
		//k = sqrt(Points.size());
		m = 2 * k + 2;
		spline();
		gauss();
		pointsSpline = getVectorPointsOutSpline();
	}
	if (pointsSpline.size() > 0) 
	{
		glLineWidth(5);
		glColor3ub(ColorRS, ColorGS, ColorBS);
		glBegin(GL_LINE_STRIP);
		for (int i = 0; i < pointsSpline.size(); i++)
		{
			glVertex2f(pointsSpline[i].x, pointsSpline[i].y);
		}
		glEnd();
		// деление на отрезки
		glLineWidth(2);
		glColor3ub(150, 50, 255);
		glEnable(GL_LINE_STIPPLE);
		glLineStipple(1, 0x00FF);
		glBegin(GL_LINES);
		glVertex2f(Points[nvtr[0][0]].x, 2 * Height);
		glVertex2f(Points[nvtr[0][0]].x, -Height);
		for (int i = 0; i < nvtr[1].size(); i++)
		{
			glVertex2f(Points[nvtr[1][i]].x, 2 * Height);
			glVertex2f(Points[nvtr[1][i]].x, - Height);
		}
		glEnd();

		glDisable(GL_LINE_STIPPLE);
	}
	// рамка с параметрами
	glColor3ub(100, 200, 200);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_DST_ALPHA);
	glBegin(GL_QUADS);
	glVertex2f(Width, 0);
	glVertex2f(Width, 5 * h);
	glVertex2f(Width - h * 7, 5 * h);
	glVertex2f(Width - h * 7, 0);
	glDisable(GL_BLEND);
	glEnd();
	// параметры
	glColor3ub(0, 0, ColorB);
	char buffer[_CVTBUFSIZE];
	// альфа
	string al = "alpha   = ";
	glRasterPos2f(Width - 6 * h, 4 * h);
	for (unsigned int j = 0; al[j] != '\0'; ++j)
	{
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, al[j]);
	}
	glRasterPos2f(Width - 3 * h , 4 * h);
	_gcvt_s(buffer, _CVTBUFSIZE, alpha, 3);

	if (buffer[2] == '\0')
		for (unsigned int j = 0; buffer[j] != ','; ++j)
		{
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, buffer[j]);
		}
	else
		for (unsigned int j = 0; buffer[j] != '\0'; ++j)
		{
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, buffer[j]);
		}
	// бета
	string b = "  beta   = ";
	glRasterPos2f(Width - 6 * h, 3 * h);

	for (unsigned int j = 0; b[j] != '\0'; ++j)
	{
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, b[j]);
	}

	glRasterPos2f(Width - 3 * h, 3 * h);
	_gcvt_s(buffer, _CVTBUFSIZE, beta, 3);
	if (buffer[2] == '\0')
		for (unsigned int j = 0; buffer[j] != ','; ++j)
		{
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, buffer[j]);
		}
	else
		for (unsigned int j = 0; buffer[j] != '\0'; ++j)
		{
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, buffer[j]);
		}
	// количество точек
	string nk = "      n   = ";
	glRasterPos2f(Width - 6 * h, 2 * h);
	for (unsigned int j = 0; nk[j] != '\0'; ++j)
	{
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, nk[j]);
	}
	glRasterPos2f(Width - 3 * h, 2 * h);
	_gcvt_s(buffer, _CVTBUFSIZE, Points.size(), 2);
	if (Points.size() < 10)
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, buffer[0]);
	else
	for (unsigned int j = 0; buffer[j] != '\0'; ++j)
	{
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, buffer[j]);
	}
	// количество конечных элементов
	string mk = "      k   = ";
	glRasterPos2f(Width - 6 * h, h);
	for (unsigned int j = 0; mk[j] != '\0'; ++j)
	{
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, mk[j]);
	}
	glRasterPos2f(Width - 3 * h, h);
	_gcvt_s(buffer, _CVTBUFSIZE, k, 2);
	if (k < 10)
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, buffer[0]);
	else
		for (unsigned int j = 0; buffer[j] != '\0'; ++j)
		{
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, buffer[j]);
		}
}

// Функция вывода на экран
void Display(void)
{
	glClearColor(1, 1, 1, 1);
	glClear(GL_COLOR_BUFFER_BIT);
	drawGrid();
	drawGridNumber();
	drawPoint();
	drawSpline();
	glFinish();
}

// Функция изменения размеров окна
void Reshape(GLint w, GLint h)
{
	for (int i = 0; i < Points.size(); i++)
	{
		Points[i].x += (w - Width) / 2;
		Points[i].y += (h - Height) / 2;
	}
	Width = w;
	Height = h;
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, w, 0, h, -1.0, 1.0);
	glMatrixMode(GL_MODELVIEW);
}

// Функция обработки сообщений от клавиатуры
void Keyboard(unsigned char key, int x, int y)
{
	// Изменение кординат графика
	if (key == 'w' && Points.size() != 0 && -poz_y < Height / 2 - h)
	{
		for (int i = 0; i < Points.size(); i++)
			Points[i].y -= 5;
		poz_y -= 5;
	}
	if (key == 's' && Points.size() != 0 && +poz_y < Height / 2 - h)
	{
		for (int i = 0; i < Points.size(); i++)
			Points[i].y += 5;
		poz_y += 5;
	}
	if (key == 'a' && Points.size() != 0 && poz_x < Width / 2 - 2 * h)
	{
		for (int i = 0; i < Points.size(); i++)
			Points[i].x += 5;
		poz_x += 5;
	}
	if (key == 'd' && Points.size() != 0 && -poz_x < Width / 2 - h)
	{
		for (int i = 0; i < Points.size(); i++)
			Points[i].x -= 5;
		poz_x -= 5;
	}
	// масштаб
	if (key == 'u' && scale < 8)
	{
		scale *= 2;
		// траектория
		for (int i = 0; i < Points.size(); i++)
		{
			// положительные координаты x
			if (Points[i].x > Width / 2 + poz_x)
			{
				Points[i].x -= Width / 2 + poz_x;
				Points[i].x /= 2;
				Points[i].x += Width / 2 + poz_x;
			}
			else // отрицательные координаты x
			{
				Points[i].x = Width / 2 + poz_x - Points[i].x;
				Points[i].x /= 2;
				Points[i].x = Width / 2 + poz_x - Points[i].x;
			}
			// положительные координаты y
			if (Points[i].y > Height / 2 + poz_y)
			{
				Points[i].y -= Height / 2 + poz_y;
				Points[i].y /= 2;
				Points[i].y += Height / 2 + poz_y;
			}
			else // отрицательные координаты y
			{
				Points[i].y = Height / 2 + poz_y - Points[i].y;
				Points[i].y /= 2;
				Points[i].y = Height / 2 + poz_y - Points[i].y;
			}
		}
	}
	
	if (key == 'i' && scale > 0.125)
	{
		scale /= 2;
		for (int i = 0; i < Points.size(); i++)
		{
			// положительные координаты x
			if (Points[i].x > Width / 2 + poz_x)
			{
				Points[i].x -= Width / 2 + poz_x;
				Points[i].x *= 2;
				Points[i].x += Width / 2 + poz_x;
			}
			else // отрицательные координаты x
			{
				Points[i].x = Width / 2 + poz_x - Points[i].x;
				Points[i].x *= 2;
				Points[i].x = Width / 2 + poz_x - Points[i].x;
			}
			// положительные координаты y
			if (Points[i].y > Height / 2 + poz_y)
			{
				Points[i].y -= Height / 2 + poz_y;
				Points[i].y *= 2;
				Points[i].y += Height / 2 + poz_y;
			}
			else // отрицательные координаты y
			{
				Points[i].y = Height / 2 + poz_y - Points[i].y;
				Points[i].y *= 2;
				Points[i].y = Height / 2 + poz_y - Points[i].y;
			}
		}
	}
	// Изменение RGB-компонент цвета сплайна (клавиатура)
	if (key == 'r')
	{
		ColorRS += 1;
	}

	if (key == 'g')
	{
		ColorGS += 1;
	}

	if (key == 'b')
	{
		ColorBS += 1;
	}
	// Изменение RGB-компонент цвета точек У ГРУППЫ (МЕНЮ)
	if (key == 'R')
	{
		ColorRS += 50;
	}

	if (key == 'G')
	{
		ColorGS += 50;
	}

	if (key == 'B')
	{
		ColorBS += 50;
	}
	// количество КЭ
	if (key == 'k' && k < Points.size() - 1)
	{
		k++;
	}

	if (key == 'l' && k > 1)
	{
		k--;
	}

	// удаление сплайна
	if (key == ' ')
	{
		Points.clear();
		pointsSpline.clear();
	}
	// изменение альфа
	if (key == '1' && alpha < 5000)
	{
		alpha += 0.5;
	}

	if (key == '2' && alpha > 0.05)
	{
		alpha -= 0.05;
	}

	// изменение бета
	if (key == '3' && beta < 5)
	{
		beta += 0.05;
	}

	if (key == '4' && beta > 0.05)
	{
		beta -= 0.05;
	}


	glutPostRedisplay();
}

// функция обработки специальных клавиш
void processSpecialKeys(int key, int x, int y)
{
	// Перемещение сплайна
	if (key == GLUT_KEY_UP)
	{
		for (int i = 0; i < Points.size(); i++)
			Points[i].y += 5;
	}

	if (key == GLUT_KEY_DOWN)
	{
		for (int i = 0; i < Points.size(); i++)
			Points[i].y -= 5;
	}

	if (key == GLUT_KEY_LEFT)
	{
		for (int i = 0; i < Points.size(); i++)
			Points[i].x -= 5;
	}

	if (key == GLUT_KEY_RIGHT)
	{
		for (int i = 0; i < Points.size(); i++)
			Points[i].x += 5;
	}

	glutPostRedisplay();
}


// Функция обработки сообщения от мыши
void Mouse(int button, int state, int x, int y)
{
	// клавиша была нажата, но не отпущена
	if (state != GLUT_DOWN) return;

	// новая точка по левому клику 
	if (button == GLUT_LEFT_BUTTON)
	{
		tmpPoint.x = x;
		tmpPoint.y = Height - y;
		Points.push_back(tmpPoint);
		// сортировка точек по значениям x
		Point tmp;
		for (int i = 0; i < Points.size() - 1; i++)
		{
			for (int j = 0; j < Points.size() - i - 1; j++)
			{
				if (Points[j].x > Points[j + 1].x) 
				{
					// меняем элементы местами
					tmp = Points[j];
					Points[j] = Points[j + 1];
					Points[j + 1] = tmp;
				}
			}
		}
	}
	glutPostRedisplay();
}

void Menu(int pos)
{
	int key = (keys)pos;

	switch (key)
	{
	case KeyW: Keyboard('w', 0, 0); break;
	case KeyS: Keyboard('s', 0, 0); break;
	case KeyA: Keyboard('a', 0, 0); break;
	case KeyD: Keyboard('d', 0, 0); break;
	case KeyU: Keyboard('u', 0, 0); break;
	case KeyI: Keyboard('i', 0, 0); break;
	case KeyR: Keyboard('R', 0, 0); break;
	case KeyG: Keyboard('G', 0, 0); break;
	case KeyB: Keyboard('B', 0, 0); break;
	case KeyS1: processSpecialKeys(GLUT_KEY_UP, 0, 0); break;
	case KeyS2: processSpecialKeys(GLUT_KEY_DOWN, 0, 0); break;
	case KeyS3: processSpecialKeys(GLUT_KEY_LEFT, 0, 0); break;
	case KeyS4: processSpecialKeys(GLUT_KEY_RIGHT, 0, 0); break;
	case Key1: Keyboard('1', 0, 0); break;
	case Key2: Keyboard('2', 0, 0); break;
	case Key3: Keyboard('3', 0, 0); break;
	case Key4: Keyboard('4', 0, 0); break;
	case KeyPr: Keyboard(' ', 0, 0); break;

	default:
		// перемещение
		int menu_move = glutCreateMenu(Menu);
		glutAddMenuEntry("Вверх", KeyW);
		glutAddMenuEntry("Вниз", KeyS);
		glutAddMenuEntry("Bлево", KeyA);
		glutAddMenuEntry("Вправо", KeyD);

		// масштаб
		int menu_scale = glutCreateMenu(Menu);
		glutAddMenuEntry("Отдалить", KeyU);
		glutAddMenuEntry("Приблизить", KeyI);

		// изменение цвета
		int menu_color = glutCreateMenu(Menu);
		glutAddMenuEntry("Компонента R", KeyR);
		glutAddMenuEntry("Компонента G", KeyG);
		glutAddMenuEntry("Компонента B", KeyB);

		// перемещение сплайна
		int menu_run = glutCreateMenu(Menu);
		glutAddMenuEntry("Вверх", KeyS1);
		glutAddMenuEntry("Вниз", KeyS2);
		glutAddMenuEntry("Bлево", KeyS3);
		glutAddMenuEntry("Вправо", KeyS4);

		// альфа
		int menu_alpha = glutCreateMenu(Menu);
		glutAddMenuEntry("Увеличить", Key1);
		glutAddMenuEntry("Уменьшить", Key2);

		// бета
		int menu_beta = glutCreateMenu(Menu);
		glutAddMenuEntry("Увеличить", Key3);
		glutAddMenuEntry("Уменьшить", Key4);

		// глобальное
		int menu = glutCreateMenu(Menu);
		glutAddSubMenu("Изменение масштаба", menu_scale);
		glutAddSubMenu("Перемещение по плоскости", menu_move);
		glutAddSubMenu("Перемещение сплайна", menu_run);
		glutAddSubMenu("Изменение цвета сплайна", menu_color);
		glutAddSubMenu("Изменение коэффициента альфа", menu_alpha);
		glutAddSubMenu("Изменение коэффициента бета", menu_beta);
		glutAddMenuEntry("Удалить сплайн", KeyPr);


		glutAttachMenu(GLUT_RIGHT_BUTTON);
		Keyboard(Empty, 0, 0);
	}
}

// Головная программа
void main(int argc, char* argv[])
{
	setlocale(LC_ALL, "Russian");
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB);
	glutInitWindowSize(Width, Height);
	glutCreateWindow("ИС вуза");
	cout << "Инструкция по использованию клавиш: " << endl << endl;

	cout << "Левая кнопка мыши - добавить новую точку" << endl;
	cout << "Пробел - удалить сплайн" << endl;
	cout << "Правая кнопка мыши - открыть меню" << endl << endl;

	cout << "u - отдаление" << endl;
	cout << "i - приближение" << endl << endl;

	cout << "k - увеличить количество КЭ" << endl;
	cout << "l - уменьшить количество КЭ" << endl << endl;

	cout << "1 - увеличить альфа" << endl;
	cout << "2 - уменьшить альфа" << endl << endl;

	cout << "3 - увеличить бета" << endl;
	cout << "4 - уменьшить бета" << endl << endl;

	cout << "r - изменение красной компоненты цвета сплайна" << endl;
	cout << "b - изменение синей компоненты цвета сплайна" << endl;
	cout << "g - изменение зеленой компоненты цвета сплайна" << endl << endl;

	cout << "wasd - перемещение по плоскости" << endl;
	cout << "0987 - перемещение сплайна" << endl;

	Menu(Empty);
	glutDisplayFunc(Display);
	glutReshapeFunc(Reshape);
	glutKeyboardFunc(Keyboard);
	glutSpecialFunc(processSpecialKeys);
	glutMouseFunc(Mouse);
	glutMainLoop();
}
