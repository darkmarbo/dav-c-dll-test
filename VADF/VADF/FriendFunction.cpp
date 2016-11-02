/*********************************************************************/
/** Name: fufill the additional functions defined in class VAD */
/*********************************************************************/
#include "VAD.h"
/* linspace function */
void linspace(RowVectorXd &Axis,double MIN_VALUE, double MAX_VALUE)
{
	if ( Axis.size() < 0 )
		throw runtime_error( "the number of the element in the axis is wrong" );
	else if ( MIN_VALUE >= MAX_VALUE )
		throw runtime_error( "the scope of the axis is wrong" );
	else
		for ( int i = 0; i < Axis.size(); i++)
			Axis[i] = MIN_VALUE+(MAX_VALUE-MIN_VALUE)*(double(i)/double(Axis.size()-1.0));
}
/* linspace_cos function */
void linspace_cos(RowVectorXd &Axis, double MIN_VALUE, double MAX_VALUE)
{
	if ( Axis.size() < 0 )
		throw runtime_error( "the number of the element in the axis is wrong" );
	else if ( MIN_VALUE >= MAX_VALUE )
		throw runtime_error( "the scope of the axis is wrong" );
	else
		for ( int i = 0; i < Axis.size(); i++ )
		{
			Axis[i] = MIN_VALUE + 2 * i;
			Axis[i] = 0.54 + 0.46 * cos(Axis[i] * PI / MAX_VALUE );
		}
}
/* fix function */
int fix( double r )
{
	if ( r < 0 ){
		double t = -r;
		return(-floor( t ) );
	}else
		return(floor( r ) );
}
/* rfft d=2,to do FFT for every row----d=2  */ 
void rfft(MatrixXd& x,int n,MatrixXcd& y)   
{
	int y_col = 1 + fix(n/2.0);
	y.resize(x.rows(),y_col);
	FFT<double> fft;
	fft.SetFlag(fft.HalfSpectrum);
	for(int ii=0;ii<x.rows();ii++)
	{
	   MatrixXd  x0   = MatrixXd::Zero(1,n);
	   x0.block(0,0,1,x.cols()) = x.row(ii);                 //add zeros operation
	   MatrixXcd temp(1,y_col);
	   fft.fwd(temp.row(0),x0.row(0));
	   y.row(ii) = temp;
	}
}
/* do irfft. every col d = 1  */
void irfft(MatrixXcd& y,int n,MatrixXd& x)   
{
	int s[2]={y.rows(),y.cols()};
	int ps  = s[0]*s[1];
	int ns  = 2;
	int d   = 1; 
	int m   = y.rows();
	int k   = y.cols();
	MatrixXcd v_temp = y;
	MatrixXcd v;
	int mm  = 1 + fix(n/2.0);
	if(mm > m)
	{
		v.resize(v_temp.rows()+(mm-m),y.cols());
		v << v_temp, MatrixXcd::Zero(mm-m,y.cols());
	}
	else if(mm < m) 
	{
		v.resize(mm,y.cols());
	    v = v_temp.block(0,0,mm,y.cols()); 
	}
	else
	{
		v = v_temp;
	}
	v_temp.resize(0,0);
	m = mm;
	v.row(m-1).imag().fill(0);
	RowVectorXd w = RowVectorXd::Ones(k);
	VectorXcd t(m);
	for ( int ii = 0; ii < m; ii++ )
	{
		t.segment(ii,1).real().fill(0.5*sin(2*PI/n*ii));
		t.segment(ii,1).imag().fill(-0.5*cos(2*PI/n*ii));
	}
	MatrixXcd z(v.rows(),v.cols());
	MatrixXcd v_flipud(v.rows(),v.cols());
	for(int ii=0;ii<v.rows();ii++)
		v_flipud.row(ii) = v.row(v.rows()-ii-1); 
	for(int jj=0;jj<v.cols();jj++)
		z.col(jj) = t.transpose().array() + 0.5;
	z = z.array()*(v_flipud.conjugate().array()-v.array()) + v.array();
	v_flipud.resize(0,0); v.resize(0,0);
	
	MatrixXcd z_temp = z; z.resize(m-1,z_temp.cols());
	z = z_temp.block(0,0,m-1,z_temp.cols());
	z_temp.resize(0,0);

    MatrixXcd zz(z.rows(),z.cols());
	FFT<double> ifft;  
	for(int ii=0;ii<z.cols();ii++)   //每一列做fft，还要补0
	{
	   MatrixXcd   x0 = MatrixXcd::Zero(1,z.rows());
	   x0.block(0,0,1,z.rows()) = z.col(ii).transpose();                 //add zeros operation
	   MatrixXcd  temp(1,z.rows()); 
	   ifft.inv(temp.row(0),x0.row(0));
	   zz.col(ii) = temp.transpose();
	}
	z.resize(0,0);
	//MatrixXd x = MatrixXd::Zero(n,k);
	x = MatrixXd::Zero(n,k);
	for(int ii=1;ii<= zz.rows();ii++)
	{
		x.row(2*(ii-1)+0) = zz.row(ii-1).real();
	    x.row(2*(ii-1)+1) = zz.row(ii-1).imag();
	}
}
/* up sorting have position of every double elements*/
void sort_up(RowVectorXd &ss,RowVectorXd &Record_sort_row,RowVectorXi &Record_sort_num_row)  
{
	double	key; 
	int     key_num ,i;
	for ( int i = 0; i <ss.size(); i++ )
	{
		Record_sort_num_row[i]	= i;
		Record_sort_row[i]	= ss[i];
	}
	for ( int j = 1;j < ss.size(); j++ )
	{
		key	= Record_sort_row[j];
		key_num = j;
		i	= j - 1;
		while ( i >= 0 && Record_sort_row[i] > key )
		{
			Record_sort_row[i + 1]		= Record_sort_row[i];
			Record_sort_num_row[i + 1]	= Record_sort_num_row[i];
			--i;
		}
		Record_sort_row[i + 1]		= key;
		Record_sort_num_row[i + 1]	= key_num;
	}
}
/* up sorting without position of every double elements*/
void sort_up_np(VectorXd ss,VectorXd &Record_sort_row)  
{
	double	key; 
	int     i;
	Record_sort_row.resize(ss.size());
	Record_sort_row	= ss;
	for ( int j = 1;j < ss.size(); j++ )
	{
		key	= Record_sort_row[j];
		i	= j - 1;
		while (i >= 0 && Record_sort_row[i] > key )
		{
			Record_sort_row[i + 1] = Record_sort_row[i];
			--i;
		}
		Record_sort_row[i + 1]	= key;
	}
}
/* up sorting without position of every int elements*/
void sort_int_up(RowVectorXi &ss,RowVectorXi &Record_sort_row)
{
	int	i;
	int	key, key_num;
	Record_sort_row = ss;
	for ( int j = 1; j < ss.size(); j++ )
	{
		key	= Record_sort_row[j];
		key_num = j;
		i	= j - 1;
		while ( i >= 0 && Record_sort_row[i] > key )
		{
			Record_sort_row[i + 1]	= Record_sort_row[i];
			i			= i - 1;
		}
		Record_sort_row[i + 1] = key;
	}
}
/* descend sorting have position of every int elements*/
void Sort_descend(VectorXd& ss,VectorXd& Record_sort_row,VectorXi &Record_sort_num_row)  
{
	int	  i, key_num;
	int   num = ss.size();  
	double	key;
	Record_sort_row = ss;
	for ( int i = 0; i < num; i++ )
		Record_sort_num_row[i]	= i;
	for ( int j = 1; j < num; j++ )
	{
		key	= Record_sort_row[j];
		key_num = j;
		i	= j - 1;
		while ( i >= 0 && Record_sort_row[i] < key )
		{
			Record_sort_row[i + 1]		= Record_sort_row[i];
			Record_sort_num_row[i + 1]	= Record_sort_num_row[i];
			i				= i - 1;
		}
		Record_sort_row[i + 1]		= key;
		Record_sort_num_row[i + 1]	= key_num;
	}
}
/* find function such as in matlab*/
/* rules:  1  ------ >      0 ------ =   -1 ------ <
		   10  ----- >=   -10 ------ <=                   */
void find(RowVectorXd& vec,double limit,int logic,RowVectorXi& pos)
{

	if(logic == 1)              // > 
	{
		RowVectorXi v_temp = (vec.array() > limit).select(RowVectorXi::Ones(vec.size()),
		       RowVectorXi::Zero(vec.size()));
		int num1 = v_temp.sum(); v_temp.resize(0);
		pos.resize(num1);
		int index = 0;
		for(int ii=0;ii<vec.size();ii++){
	      if(vec(ii) > limit)  
			  pos(index++) = ii+1;
		}
	}
	else if(logic == 10)       // >=
	{
		RowVectorXi v_temp = (vec.array() >= limit).select(RowVectorXi::Ones(vec.size()),
		       RowVectorXi::Zero(vec.size()));
		int num1 = v_temp.sum(); v_temp.resize(0);
		pos.resize(num1);
		int index = 0;
		for(int ii=0;ii<vec.size();ii++){
	      if(vec(ii) >= limit)  
			  pos(index++) = ii+1;
		}
	}

	else if(logic == -1)       // <
	{
		RowVectorXi v_temp = (vec.array() < limit).select(RowVectorXi::Ones(vec.size()),
		       RowVectorXi::Zero(vec.size()));
		int num1 = v_temp.sum(); v_temp.resize(0);
		pos.resize(num1);
		int index = 0;
		for(int ii=0;ii<vec.size();ii++){
	      if(vec(ii) < limit)  
			  pos(index++) = ii+1;
		}
	}
	else if(logic == -10)       //  <=
	{
		RowVectorXi v_temp = (vec.array() <= limit).select(RowVectorXi::Ones(vec.size()),
		       RowVectorXi::Zero(vec.size()));
		int num1 = v_temp.sum(); v_temp.resize(0);
		pos.resize(num1);
		int index = 0;
		for(int ii=0;ii<vec.size();ii++){
	      if(vec(ii) <= limit)  
			  pos(index++) = ii+1;
		}
	}
	else if(logic == 0)       //  <=
	{
		RowVectorXi v_temp = (vec.array() == limit).select(RowVectorXi::Ones(vec.size()),
		       RowVectorXi::Zero(vec.size()));
		int num1 = v_temp.sum(); v_temp.resize(0);
		pos.resize(num1);
		int index = 0;
		for(int ii=0;ii<vec.size();ii++){
	      if(vec(ii) <= limit)  
			  pos(index++) = ii+1;
		}
	}
	else
        return;
}
/* find function with two conditions and & */
void find_and(RowVectorXd& vec,double low,double high,RowVectorXi& pos)
{
	RowVectorXi v_temp = (vec.array()>low && vec.array()<high).
		select(RowVectorXi::Ones(vec.size()),RowVectorXi::Zero(vec.size()));
	int num1 = v_temp.sum(); v_temp.resize(0);
	pos.resize(num1);
	int index = 0;
	for(int ii=0;ii<vec.size();ii++){
	      if((vec(ii)>low) && (vec(ii)<high))  
			  pos(index++) = ii+1;
		}
}
/*Get the eigen value and vector of a 2*2 Matrix */
void eig(Matrix2d& v_,Matrix2d& uvk_,Matrix2d& dvk_)	
{
	double	x1	= v_(0,0), x2 = v_(0,1);
	double	y1	= v_(1,0), y2 = v_(1,1);

	double	lamta1	= ( (x1 + y2) - sqrt( (x1 + y2) * (x1 + y2) - 4 * (x1 * y2 - x2 * y1) ) ) / 2.0;
	double	lamta2	= ( (x1 + y2) + sqrt( (x1 + y2) * (x1 + y2) - 4 * (x1 * y2 - x2 * y1) ) ) / 2.0;
	dvk_(0,0) = lamta1; dvk_(0,1) = 0;
	dvk_(1,0) = 0;      dvk_(1,1) = lamta2;
	double	T1_1	= 1;
	double	T2_1	= (lamta1 - x1) / x2;
	double	v1_1	= -1 * T1_1 / sqrt( T1_1 * T1_1 + T2_1 * T2_1 );
	double	v2_1	= -1 * T2_1 / sqrt( T1_1 * T1_1 + T2_1 * T2_1 );
	uvk_(0,0) = v1_1;
	uvk_(1,0) = v2_1;

	double	T1_2	= 1;
	double	T2_2	= (lamta2 - x1) / x2;
	double	v1_2	= T1_2 / sqrt( T1_2 * T1_2 + T2_2 * T2_2 );
	double	v2_2	= T2_2 / sqrt( T1_2 * T1_2 + T2_2 * T2_2 );
	uvk_(0,1) = v1_2;
	uvk_(1,1) = v2_2;
}
/* get the median value of a array*/
double median(VectorXd &Array)
{
	/* insert sort*/
	if (Array.size() == 0)
		return NULL;
	VectorXd Array_temp = Array;
	int    i,j;
	double temp; 
	double m_value;
	for(i = 1;i<Array_temp.size();i++)
	{
		temp = Array_temp[i];
		for(j=i-1; j>0 && temp<Array_temp[j];j--){
			Array_temp[j+1] = Array_temp[j];
		}
		Array_temp[j+1] = temp;
	}
	if(Array_temp.size() % 2 == 0)
		m_value = (Array_temp[Array_temp.size()/2-1] + Array_temp[Array_temp.size()/2])/2;
	else
	    m_value = Array_temp[Array_temp.size()/2];
	return m_value; 
}
/* hanning window */
void hamming( int N, RowVectorXd& win)
{
	for ( int ii = 0; ii < (N + 1) / 2; ii++ )
	{
		win[ii]		= 0.54 - 0.46 * cos( 2 * PI * ii / (N - 1.0) );
		win[N - 1 - ii] = win[ii];
	}
}