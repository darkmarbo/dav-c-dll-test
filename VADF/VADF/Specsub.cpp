/*-------------------------------------------------------------------*/
/**  Name: fufill specsub rats noise_seg model	                    **/
/*-------------------------------------------------------------------*/
#include "VAD.h"

void VAD::specsub_rats_noiseseg_1fn(VectorXd& si,VectorXd& noise_seg,VectorXd &pv01,VectorXd& ss_end)
{
	double		fs	= _fs;
	Sub_Algo_param	qq;             /* default algorithm constants */

	qq.of	= 2;
	qq.ti	= 16e-3;
	qq.ri	= 0;
	qq.g	= 1;
	qq.e	= 1;
	qq.am	= 3;
	qq.b	= 0.01;
	qq.al	= -5;
	qq.ah	= 20;
	qq.bt	= -1;
	qq.ne	= 0;
	qq.mx	= 0;
	qq.gh	= 1;
	qq.tf	= 'g';
	qq.rf	= 0;

	int	    ni	    = int(qq.ti * fs + 0.5);             /* frame increment in samples */
	double	tinc	= ni / fs;
	double	ne	    = qq.ne;
	int	    no		= int(qq.of + 0.5);
	int	    nf		= ni * no; 

	RowVectorXd w1 = RowVectorXd::Zero(nf+1);
	hamming(w1.size(),w1);
    RowVectorXd w  = w1.head(w1.size()-1).array().sqrt();  w1.resize(0);
	int	 N_sum	= 1 + (nf - 1) / ni;
	RowVectorXd w_sum(N_sum);
	for( int ii = 0;ii<N_sum; ii++)   w_sum(ii) = w(0+ii*ni);
	w = w.array()/ sqrt( w_sum.array().square().sum());
	w_sum.resize(0);
	MatrixXd y; VectorXd tt;
	enframe_v2(si,w,ni,y,tt);
	MatrixXcd yf;        
	rfft(y,nf,yf);
	MatrixXcd yf_temp = yf;
	yf = yf.array() * yf.conjugate().array();
    MatrixXd yp;  
	yp = yf.real();  yf.resize(0,0);
	int nr  = yp.rows();
	int nf2 = yp.cols();
	if(noise_seg.size() < nr)
	{
		VectorXd noise_seg_temp = noise_seg;
		noise_seg.resize(nr);
		noise_seg << noise_seg_temp,0;
	}
	MatrixXd dp;
	estnoisem_rats_noiseseg(yp,noise_seg,pv01,tinc,dp);
	VectorXd ssv = VectorXd::Zero(ni*(no-1));
	// without if ~nr  
	MatrixXi mz  = (yp.array() == 0).select(MatrixXi::Ones(yp.rows(),yp.cols()),MatrixXi::Zero(yp.rows(),yp.cols()));
	VectorXd ypf = yp.rowwise().sum();
	VectorXd dpf = dp.rowwise().sum();
	VectorXi mzf = (dpf.array() == 0).select(VectorXi::Ones(dpf.size()),VectorXi::Zero(dpf.size()));
	VectorXd af_temp = ypf.array()/(dpf.cast<double>().array()+mzf.cast<double>().array());
	for(int ii=0;ii<af_temp.size();ii++)
	{	
		af_temp(ii)=1+(qq.am-1)*(min(max(10*log10(af_temp(ii)),qq.al),qq.ah)-qq.ah)/(qq.al-qq.ah);
		if(mzf(ii) == 1){  af_temp(ii) = 1;}
	}
	MatrixXd v  = (dp.array()/(yp.array()+mz.cast<double>().array())).sqrt();
	dp.resize(0,0); mz.resize(0,0);
	af_temp = af_temp.array().sqrt();
	double bf = sqrt(qq.b);
	MatrixXd  af(af_temp.size(),nf2);
	for(int ii=0;ii<nf2;ii++)
		af.col(ii) = af_temp;
	af_temp.resize(0);
	MatrixXi mf = (v.array() >= (af.array()+bf).inverse()).select
                  (MatrixXi::Ones(af.rows(),af.cols()),MatrixXi::Zero(af.rows(),af.cols()));
	MatrixXd g  = MatrixXd::Zero(v.rows(),v.cols());
	double  eg  = qq.e/qq.g;
	double  gh  = qq.gh;
	for(int ii=0;ii<mf.rows();ii++)
	{
		for(int jj=0;jj<mf.cols();jj++)
		{
			if(mf(ii,jj) == 1){ g(ii,jj) = min(bf*v(ii,jj),gh); }
			else              { g(ii,jj) = 1-af(ii,jj)*v(ii,jj);}
		}
	}
	g = qq.mx + (1-qq.mx)*g.array();
	MatrixXcd out(yf_temp.rows(),yf_temp.cols());
	out.real() = yf_temp.real().array()* g.array();
	out.imag() = yf_temp.imag().array()* g.array() ; yf_temp.resize(0,0);
	MatrixXcd out_p_temp = out.array() * out.conjugate().array();
	MatrixXd  out_p      = out_p_temp.real(); out_p_temp.resize(0,0);
	VectorXd  out_pf     = out_p.rowwise().sum();

	for(int i=1;i<=nr;i++)
	{
		if(yp.block(i-1,0,1,7).sum() > (yp.block(i-1,0,1,yp.cols()).sum())/2.0)
		{
			yp.block(i-1,0,1,7).fill(0);
		    out_p.block(i-1,0,1,7).fill(0);
			out.block(i-1,0,1,7).fill((0,0));	
		}
	}
	VectorXd out_dpft  = 0.5*dpf;
	VectorXi out_smth  = VectorXi::Ones(nr);
	for(int i=1;i<=nr;i++)
	{
		if(out_pf(i-1) < out_dpft(i-1)) out_smth(i-1) = 0;
		else out_smth(i-1) = 1;
	}
	for(int i=1;i<=10;i++)
	{
		if(out_pf(i-1) < out_dpft(i-1)) 
			out.row(i-1).fill((0,0)); 
	}
	for(int i=nr-9;i<= nr;i++)
	{
		if(out_pf(i-1) < out_dpft(i-1)) 
			out.row(i-1).fill((0,0)); 
	}
	for(int i=11;i<=nr-10;i++)
	{
		if((out_pf(i-1) < out_dpft(i-1))&&(out_smth.segment(i-11,21).sum()<5))   //;,2015.8.24.23:18 ÕÒµ½Õâ¸öbitch ;  
			out.row(i-1).fill((0,0));
	}
	MatrixXd se_temp;						  //save irfft output
	MatrixXcd out_temp = out.transpose();       //need to do irfft 	
	irfft(out_temp,nf,se_temp);
	out_temp.resize(0,0); out.resize(0,0);
	MatrixXd se = se_temp.transpose();   se_temp.resize(0,0);
	for(int ii=0;ii<se.rows();ii++)
	{
		se.row(ii) = se.row(ii).array() * w.array(); 
	}
	int n1,n2;   //n1:0,2,4,6  n2:1,3,5,7
	if(nr % 2 == 0){n1 = n2 = nr/2;}
	else{ n1 = nr/2+1; n2 =  nr/2;}
	VectorXd se_1(n1*se.cols()); VectorXd se_2(n2*se.cols());
	for(int ii=0;ii<nr;ii++)
	{
		if(ii%2 == 0)
			se_1.segment(ii/2*se.cols(),se.cols()) = se.row(ii);
		else
			se_2.segment(ii/2*se.cols(),se.cols()) = se.row(ii);
	}
	se.resize(0,0);
	MatrixXd ss = MatrixXd::Zero(ni*(nr+no-1),no);
	for(int ii=1;ii<=no;ii++)
	{
		int nm = nf * (1 + floor(double(nr-ii)/no));
	    if(ii==1){ss.block( (ii-1)*ni,ii-1,nm,1) = se_1;}
		else {ss.block( (ii-1)*ni,ii-1,nm,1) = se_2;}
	}
	VectorXd ss_sum = ss.rowwise().sum();
	ss_end = ss_sum.segment(0,si.size());  //the end output
	//cout << ss_end.segment(736,4) << endl;
}

/* may be exist some bug,but it was only noise estnoisem,look t's output*/
void VAD::estnoisem_rats_noiseseg(MatrixXd &yf,VectorXd& noise_seg,VectorXd &pv01,double tz,MatrixXd &x)
{
	int nr      = yf.rows(); 
	int nrf     = yf.cols();
	VectorXd strongN_n  = VectorXd::Zero(nr);
	VectorXd strongN_nn = VectorXd::Zero(nr);
	VectorXd yf_strongN(yf.rows());	
	yf_strongN = yf.block(0,79,yf.rows(),50).rowwise().sum();
	for(int i=1;i<= nr;i++)
	{
		if(yf_strongN(i-1) > 5) 
			strongN_n(i-1) = 1;
	}
	for(int i=1;i<=nr-10;i++)
	{
		if ( strongN_n.segment(i-1,11).sum() > 6)
			strongN_nn.segment(i-1,11).fill(1);
	}
	//MatrixXd x  = MatrixXd::Zero(nr,nrf); 
	x  = MatrixXd::Zero(nr,nrf); 
	double		tinc	= tz;
	int		    nrcum	= 0;
	Est_Algo_param	qq;
	qq.taca		= 0.0449;                                                               /* smoothing time constant for alpha_c = -tinc/log(0.7) in equ (11) */
	qq.tamax	= 0.392;                                                                /* max smoothing time constant in (3) = -tinc/log(0.96) */
	qq.taminh	= 0.0133;                                                               /* min smoothing time constant (upper limit) in (3) = -tinc/log(0.3) */
	qq.tpfall	= 0.064;                                                                /* time constant for P to fall (12) */
	qq.tbmax	= 0.0717;                                                               /* max smoothing time constant in (20) = -tinc/log(0.8) */
	qq.qeqmin	= 2;                                                                    /* minimum value of Qeq (23) */
	qq.qeqmax	= 14;                                                                   /* max value of Qeq per frame */
	qq.av		= 2.12;                                                                 /* fudge factor for bc calculation (23 + 13 lines) */
	qq.td		= 1.536;                                                                /* time to take minimum over */
	qq.nu		= 8;                                                                    /* number of subwindows */
	qq.qith[0]	= 0.03; qq.qith[1] = 0.05;
	qq.qith[2]	= 0.06; qq.qith[3] = MAX_DBL;                                           /* noise slope thresholds in dB/s */
	qq.nsmdb[0]	= 47;   qq.nsmdb[1] = 31.4;
	qq.nsmdb[2]	= 15.7; qq.nsmdb[3] = 4.1;

	/* unpack paramter structure */
	double	taca		= qq.taca;                                                      /* smoothing time constant for alpha_c = -tinc/log(0.7) in equ (11) */
	double	tamax		= qq.tamax;                                                     /* max smoothing time constant in (3) = -tinc/log(0.96) */
	double	taminh		= qq.taminh;                                                    /* min smoothing time constant (upper limit) in (3) = -tinc/log(0.3) */
	double	tpfall		= qq.tpfall;                                                    /* time constant for P to fall (12) */
	double	tbmax		= qq.tbmax;                                                     /* max smoothing time constant in (20) = -tinc/log(0.8) */
	double	qeqmin		= qq.qeqmin;                                                    /* minimum value of Qeq (23) */
	double	qeqmax		= qq.qeqmax;                                                    /* max value of Qeq per frame */
	double	av			= qq.av;                                                        /* fudge factor for bc calculation (23 + 13 lines) */
	double	td			= qq.td;                                                        /* time to take minimum over */
	double	nu			= qq.nu;                                                        /* number of subwindows */
	double	qith[4]		= { qq.qith[0], qq.qith[1], qq.qith[2], qq.qith[3] };           /* noise slope thresholds in dB/s */
	double	nsmdb[4]	= { qq.nsmdb[0], qq.nsmdb[1], qq.nsmdb[2], qq.nsmdb[3] };;      /* maximum permitted +ve noise slope in dB/s */

	double	aca			= exp( -1 * tinc / taca );
	double	acmax		= aca;                                                                  /* min value of alpha_c = 0.7 in equ (11) also = 0.7 */
	double	amax		= exp( -1 * tinc / tamax );                                             /* max smoothing constant in (3) = 0.96 */
	double	aminh		= exp( -1 * tinc / taminh );                                            /* min smoothing constant (upper limit) in (3) = 0.3 */
	double	bmax		= exp( -1 * tinc / tbmax );                                             /* max smoothing constant in (20) = 0.8 */
	double	snrexp		= -1 * tinc / tpfall;
	int	nv	= int(td /(tinc * nu) + 0.5);                                                  /* length of each subwindow in frames */
	if (nv < 4)
	{
		nv	= 4; nu	= max(int(td/(tinc*nv)+0.5),1);
	}
	double	md, hd;
	double	mv, hv;
	int	nd = nu * nv;
	mhvals( nd, md, hd );
	mhvals( nv, mv, hv );
	double nsms[4] = {};
	for ( int ii = 0; ii < 4; ii++ )
		nsms[ii] = pow(10,nsmdb[ii]*nv*tinc/10.0);
	double	qeqimax = 1.0 / qeqmin;
	double	qeqimin = 1.0 / qeqmax;

	RowVectorXd p     = yf.row(0);
	RowVectorXd ne_min = p;
	double      ac    = 1;  
	RowVectorXd sn2   = p;
	RowVectorXd pb    = p;
	RowVectorXd pb2   = pb.array().square();
	RowVectorXd pminu = p;
	RowVectorXd actmin(nrf); actmin.fill(MAX_DBL);
	RowVectorXd actminsub = actmin;
	int	        subwc = nv;
	MatrixXd    actbuf(int(nu),nrf); actbuf.fill(MAX_DBL);
	int         ibuf  = 0;
	RowVectorXi lminflag = RowVectorXi::Zero(nrf);

	MatrixXd    yf_temp = yf.block(0,0,min(50,nr),yf.cols());
	if( pv01.head(10).sum() >= 1)
	{
		for(int ii=1;ii<= yf_temp.cols();ii++)
			p(ii-1) = yf_temp.col(ii-1).minCoeff();
	    ne_min = p; sn2 = p; pb = p; pminu = p;
		pb2    = pb.array().square();
	}
	RowVectorXd yft(yf.cols());
	RowVectorXd ah(p.size());
	RowVectorXd b(p.size());
	RowVectorXd qeqi(pb.size());
	RowVectorXd qeqi_temp(pb.size());
	RowVectorXd bmind(qeqi.size());
	RowVectorXd bminv(qeqi.size());
	RowVectorXi kmod(p.size());
	RowVectorXi lmin(kmod.size());
	RowVectorXi i;
	for(int t=1; t<=nr;t++)
	{
		yft  = yf.row(t-1);
	    double acb = 1/(1+ (p.sum()/yft.sum()-1)*(p.sum()/yft.sum()-1));
        ac  = aca*ac+(1-aca)* max(acb,acmax);
		ah = amax * ac * ((1+ (p.array()/sn2.array()-1).square()).inverse());
		double snr = p.sum()/sn2.sum();
		for (int jj =0;jj < ah.size(); jj++ )
			ah[jj] = max(ah[jj], min(aminh,pow(snr,snrexp)));
		if((t>1 && strongN_nn(t-1) !=0)||(noise_seg(t-1) != 0)||(t<11 && pv01.head(10).sum()>=1))
		{
			;//cout << t << endl; ok
		}
		else
		{
			p = ah.array() * p.array() + (1-ah.array()) * yft.array();
		    for(int jj=0;jj<ah.size();jj++) b(jj) = min( ah(jj) * ah(jj), bmax);
			pb  = b.array()*pb.array() + (1-b.array())*p.array();
			pb2 = b.array()*pb2.array() + (1-b.array())* (p.array().square());
		}
		qeqi_temp = (pb2.array()-pb.array().square())/(2*sn2.array().square());
		for(int jj=0;jj<qeqi_temp.size();jj++)
			qeqi(jj) = max(min(qeqi_temp(jj),qeqimax),qeqimin/(t+nrcum));
		double  qiav = qeqi.sum()/nrf;
		double  bc   = 1 + av*sqrt(qiav); 
		bmind = 1+2*(nd-1)*(1-md)/(qeqi.array().inverse()-2*md);
		bminv = 1+2*(nv-1)*(1-mv)/(qeqi.array().inverse()-2*mv);
		for(int jj=0;jj<p.size();jj++)
			kmod(jj) = (bc*p(jj)*bmind(jj) < actmin(jj))? 1:0;
		if(kmod.any())
		{
			//cout << t << endl;  ok
			for (int jj=0;jj < p.size();jj++ )
			{
				if ( kmod(jj) != 0 )
				{
					actmin(jj)	  = bc * p(jj) * bmind(jj);
					actminsub(jj) = bc * p(jj) * bminv(jj);
				}
			}
		}
		if(subwc>1 && subwc<nv)
		{
			for(int jj=0;jj<lminflag.size();jj++ )
			{
				lminflag(jj) = lminflag(jj) | kmod(jj);
				pminu(jj)	 = min( actminsub(jj),pminu(jj));
				sn2(jj)		 = pminu(jj);
			}
		}
		else
		{
			if( subwc >= nv)
			{
				ibuf = 1+ibuf-(fix(ibuf/nu)*nu); 
			    actbuf.row(ibuf-1) = actmin;
			    for(int jj=0;jj<actbuf.cols();jj++)
					pminu(jj) = actbuf.col(jj).minCoeff();
			    
				int ni =0,n=0;
				for(int ii=0;ii<4;ii++ ){if(qiav<qith[ii]) ni++;}
				i.resize(ni);
				for(int ii=0;ii<4;ii++ ){n=0;if(qiav<qith[ii]){i[n]=ii+1;n++;}}
			    double nsm = nsms[i(0)-1];
			    for(int ii = 0;ii< kmod.size();ii++ )
					lmin[ii] = lminflag[ii] & (kmod[ii] == 0)&(actminsub[ii] < nsm * pminu[ii])&(actminsub[ii] > pminu[ii]);
			    if(lmin.any())
				{
					//cout << t << endl;
					for(int jj=0; jj<pminu.size(); jj++)
					{
						if ( lmin[jj] != 0 )
						{
							pminu[jj] = actminsub[jj];
							actbuf.col(jj).fill(pminu[jj]);
						}
					}
				}
				lminflag.fill(0);
				actmin.fill(MAX_DBL);
				subwc = 0; 
			}
		}
		++subwc;
		x.row(t-1)  = sn2;
	}
	int ii=0;
}
void VAD::mhvals( int d, double & m_, double & h_ )
{
	double	dmh[54] = {
		1,   0,	    0,
		2,   0.26,  0.15,
		5,   0.48,  0.48,
		8,   0.58,  0.78,
		10,  0.61,  0.98,
		15,  0.668, 1.55,
		20,  0.705, 2,
		30,  0.762, 2.3,
		40,  0.8,   2.52,
		60,  0.841, 3.1,
		80,  0.865, 3.38,
		120, 0.89,  4.15,
		140, 0.9,   4.35,
		160, 0.91,  4.25,
		180, 0.92,  3.9,
		220, 0.93,  4.1,
		260, 0.935, 4.7,
		300, 0.94,  5
	};            /* dmh[18*3] */
	int	i_num = 0;
	for ( int ii = 0; ii < 18; ii++ )
	{
		if ( d <= dmh[ii * 3] )
			i_num++;
	}
	int	* i_i	= new int[i_num]();
	int	n	= 0;
	for ( int ii = 0; ii < 18; ii++ )
	{
		if ( d <= dmh[ii * 3] )
			i_i[n++] = ii + 1;
	}
	int	i, j;
	double	m, h;
	double	qj, qi, q;
	/* get the paramter of i and j */
	if ( i_num == 0 )
	{
		i	= 18;
		j	= i;
	}else  {
		i	= i_i[0];
		j	= i - 1;
	}
	/* the next if */
	if ( d == dmh[(i - 1) * 3] )
	{
		m	= dmh[(i - 1) * 3 + 1];
		h	= dmh[(i - 1) * 3 + 2];
	}else  {
		qj	= sqrt( dmh[(i - 2) * 3] );
		qi	= sqrt( dmh[(i - 1) * 3] );
		q	= sqrt( double(d) );
		h	= dmh[(i - 1) * 3 + 2] + (q - qi) * (dmh[(j - 1) * 3 + 2] - dmh[(i - 1) * 3 + 2]) / (qj - qi);
		m	= dmh[(i - 1) * 3 + 1] + (qi * qj / q - qj) * (dmh[(j - 1) * 3 + 1] - dmh[(i - 1) * 3 + 1]) / (qi - qj);
	}
	m_	= m;
	h_	= h;
	delete[]i_i; i_i = NULL;
}

/* enframe */
void VAD::enframe_v2(VectorXd& x,RowVectorXd& win,int inc,MatrixXd &f,VectorXd &t)
{
	int nx     = x.size();
	int nwin   = win.size();
	int lw     = nwin;
	VectorXd w = win.transpose();	
	int nli    = nx-lw+inc;  
	int nf     = fix(double(nli)/inc);
	int na     = nli - inc * nf;
		f      = MatrixXd::Zero(nf+1,lw);
	VectorXi    indf =  VectorXi::Zero(nf);
	RowVectorXi inds =  RowVectorXi::Zero(lw);
	for(int ii = 0; ii < nf; ii++ )  indf[ii] = inc * ii;
	for(int ii = 0; ii < lw; ii++ )  inds[ii] = ii + 1;                             	
	for(int ii = 0; ii < nf; ii++ )
	{
		for (int jj = 0; jj < lw; jj++ )
			f(ii,jj) = w(jj) * x(inds(jj) + indf(ii) - 1);
	}
	RowVectorXi ix(lw);
	for (int ii = 0; ii<ix.size();ii++)  ix(ii)=1+(nx-na+ii)%(2*nx);
	int* nnx = new int[ix.size()]();
	for (int ii=0;ii<ix.size();ii++ )
	{
		if ( ix[ii] > nx )
			nnx[ii] = ix[ii] + 1 * (2 * nx + 1 - 2 * ix[ii]);
		else
			nnx[ii] = ix[ii] + 0 * (2 * nx + 1 - 2 * ix[ii]);
		f(nf,ii) = x[nnx[ii] - 1] * w[ii];
	}
	delete[]nnx;  nnx = NULL;
	double	t0= (1 + lw) / 2.0;
		    t = VectorXd::Zero(nf);
	for(int ii=0;ii<nf;ii++)
	    t[ii] = t0 + inc * ii;
}



