/*-------------------------------------------------------------------*/
/**  Name: fufill pitchestm	estimation model	                    **/
/*-------------------------------------------------------------------*/
#include "VAD.h"

void VAD::pitchestm(VectorXd& data,VectorXd& pv01,VectorXd& fx)
{
	VectorXd  pv;
	fxpefac(data,fx,pv);
	int npv = pv.size();
	int nstart=0,nstop=0;
	pv01 = VectorXd::Zero(_nfr10);
	int sign_pv = 0;
	for(int ii=1;ii<=npv;ii++)
	{
		if(pv(ii-1) > 0.25)
		{
			pv01(ii+2) = 1;
			if(sign_pv == 0)
			{
				sign_pv = 1;
			    nstart  = ii;
			}	
		}
		else
		{
			if (sign_pv == 1)
			{
				sign_pv = 0;
				nstop   = ii-1;
				if(nstop - nstart <3)
				{
					pv01.segment(nstart+2,nstop-nstart+1) = VectorXd::Zero(nstop-nstart+1);
				}
			}
		}
	}
	pv01.head(3)         = pv01(3)*VectorXd::Ones(3);
	RowVectorXd fxtmp(_nfr10);
	fxtmp.head(3)        = fx(3)*VectorXd::Ones(3);
    fxtmp.segment(3,npv) = fx.head(npv);
	if((npv+3)< _nfr10)
	{
		pv01.segment(npv+3,(_nfr10-npv-3))  = pv01(npv+2)*VectorXd::Ones(_nfr10-npv-3);
	    fxtmp.segment(npv+3,(_nfr10-npv-3)) = fx(npv-1)*VectorXd::Ones(_nfr10-npv-3);
	}
	else
	{
		pv01  = pv01.head(_nfr10);
	    fxtmp = fxtmp.head(_nfr10); 
	}
    fx = fxtmp;
}  //VectorXd fx ,VectorXd pv
void VAD::fxpefac(VectorXd& s,VectorXd& fx,VectorXd& pv)
{
	Vector6d w_u;
	w_u <<  0.2123723,0.207788,0.2701817,0.1293616,0.04741722,0.1328791;
	MatrixXd m_u(6,2);
	m_u <<  0.2220388, 0.4067706, 0.04567656, 0.4016914, 0.8415278, 0.3192158,
		    0.2194808, 0.1910079, 1.6347,     0.5819833, 1.181519,  0.6996485;
	MatrixXd v_u(4,6);
	v_u <<  0.01413822,  0.0009377269, 0.1233703,   0.01779449,  1.110173,   0.5477135,
		    0.003357913, 0.0006220489, 0.004299293, 0.002078821, 0.00718649, -0.00182316, 
		    0.003357913, 0.0006220489, 0.004299293, 0.002078821, 0.00718649, -0.00182316,
		    0.01786169,  0.03422057,   0.007660504, 0.001605052, 0.005734435, 0.05659796;
	Vector6d w_v;
	w_v <<  0.07758689, 0.2109879, 0.1856225, 0.06853158, 0.2701563, 0.1871148;
	MatrixXd m_v(6,2);
	m_v <<  1.208656,    0.3365564, 1.216643,  0.5971916,  4.08585,   1.240948,
		    8.322102,    1.349939,  1.734108,  1.168643,   0.5107205, 0.940308;
	MatrixXd v_v(4,6);
	v_v <<  0.06181574,  0.2946077,  2.508473,    14.17252,    0.5834894,   0.05978017,
		    0.002950501, 0.01433284, -0.03310555, -0.09009174, -0.07854027, 0.005528601,
		    0.002950501, 0.01433284, -0.03310555, -0.09009174, -0.07854027, 0.005528601,	
		    0.004528442, 0.02684239, 0.1098579,   0.07989255,  0.1108958,   0.1309329 ;
	RowVector6d dpwtdef;
	dpwtdef <<  1.0000, 0.8250, 0.01868, 0.006773, 98.9, -0.4238;  
	/* Algorithm parameter defaults */
	Algo_param p;
	p.fstep		= 5;       p.fmax		= 4000;
	p.fres		= 20;      p.fbanklo	= 40;
	p.mpsmooth	= 201;     p.shortut	= 7;
	p.pefact	= 1.5;     p.numopt	    = 3;
	p.w         = dpwtdef; p.flim << 60,400;    
	p.tmf		= 2;       p.tinc		= 0.01;
	VectorXd tx; RowVectorXd f; MatrixXd MIX;
	spgrambw(s,tx,f,MIX);
	int     nframes = tx.size();
	double  txinc   = tx(1)-tx(0);
	SparseMatrix<double> trans;  RowVectorXd cf;
    filtbankm_v2(f.size(),2*f.size()-1,2*f(f.size()-1),p.fbanklo,f(f.size()-1),trans,cf);
	MatrixXd OO = MIX * trans.adjoint();
	MIX.resize(0,0);  trans.resize(0,0);

	RowVectorXd ltass;
	stdspectrum(cf,ltass);
    RowVectorXd  auxf1(1+cf.size());
	   auxf1 << cf(0),(cf.segment(0,cf.size()-1)+cf.segment(1,cf.size()-1))/2,cf.tail(1); 
	RowVectorXd  auxf(auxf1.size()-1);
	   auxf  <<  auxf1.segment(1,auxf.size()) -  auxf1.segment(0,auxf.size());
	auxf1.resize(0);
	ltass = ltass.array() * auxf.array();
	for(int jj=0;jj<OO.rows();jj++)                       // block 不能为左值
		OO.row(jj) = OO.row(jj).array() * auxf.array();  
	auxf.resize(0);
	if(tx(tx.size()-1) < p.shortut)
	{
		RowVectorXd eltass(OO.cols());
	    for(int ii=0; ii<eltass.size(); ii++) eltass(ii) = OO.col(ii).mean();
		smooth(eltass,p.mpsmooth);
	    eltass = eltass.transpose();
		double cte = ltass.mean()/eltass.mean();
		eltass = eltass.array() * cte;
		OO     = OO.array() * cte;
		RowVectorXd alpha = ltass.array()/eltass.array();
		alpha  =  alpha.transpose();
		for(int ii=0;ii<OO.rows();ii++)
			OO.row(ii) = OO.row(ii).array() * alpha.array(); 
	}
	else
	{
		double tsmo = 2;
		int	   stt	= floor(tsmo/txinc + 0.5);
		VectorXd filttime(stt+stt-1);
		   filttime << VectorXd::Ones(stt), VectorXd::Zero(stt-1);
        MatrixXd eltass = MatrixXd::Zero(OO.rows(),OO.cols());
		imfilter(OO,filttime,eltass,1);   // for column filter		
	    VectorXd filtfreq    = VectorXd::Ones(p.mpsmooth);
		MatrixXd eltass_temp = MatrixXd::Zero(OO.rows(),OO.cols());
		imfilter(eltass,filtfreq,eltass_temp,0);
		eltass = eltass_temp;  eltass_temp.resize(0,0);
		VectorXd cte(eltass.rows());
		double lta_mean = ltass.mean();
		for(int ii=0;ii<eltass.rows();ii++)
			cte(ii) = lta_mean / eltass.row(ii).mean();
		for(int ii=0;ii<eltass.cols();ii++){
			eltass.col(ii) = eltass.col(ii).array()*cte.array(); 
		    OO.col(ii)     = OO.col(ii).array()* cte.array();
		}
		MatrixXd alpha(OO.rows(),OO.cols());
		for(int ii=0;ii<eltass.rows();ii++)
			alpha.row(ii) = ltass.array() / eltass.row(ii).array();
		OO = OO.array() * alpha.array();
	}
	ltass.resize(0);

    RowVectorXi  ini;  find(cf,2*cf(0),1,ini);
	RowVectorXd  sca_temp = cf.array()/cf(ini(0)-1);
	RowVectorXi  sca_01   = (sca_temp.array()>0.5 && sca_temp.array()<10.5).
	             select(RowVectorXi::Ones(sca_temp.size()),RowVectorXi::Zero(sca_temp.size()));
	int sca_num = sca_01.sum();
	int index   = 0;
	RowVectorXd  sca(sca_num);
	for(int ii=0; ii<sca_temp.size(); ii++){
		if( sca_temp(ii)<10.5 && sca_temp(ii)> 0.5)
			sca(index++) = sca_temp(ii);}
	sca_temp.resize(0); sca_01.resize(0);
	RowVectorXd  filh = p.pefact-(2*PI*sca.array()).cos();
	for(int ii=0;ii<filh.size();ii++)
		filh(ii) = -1*log10(filh(ii));
	filh = filh.array() -  filh.mean();
	RowVectorXi  posit;  find(sca,1.0,10,posit);
	if ( posit.size() % 2 == 0)
	{
		RowVectorXd filh_copy = filh;
		filh.resize(filh_copy.size()+1);
		filh << filh_copy,0;
	}
	RowVectorXi  negat;  find(sca,1.0,-1,negat);
	int numz = posit.size() -1 - negat.size();
	filh = filh.array()/filh.maxCoeff();
	RowVectorXd filh_copy = filh;
	filh.resize(numz+filh_copy.size());
		filh << RowVectorXd::Zero(numz),filh_copy;
	filh_copy.resize(0);  negat.resize(0); sca.resize(0); ini.resize(0); posit.resize(0);
    VectorXd filh2 = filh.adjoint(); 
	
	MatrixXd B = MatrixXd::Zero(OO.rows(),OO.cols());
	imfilter(OO,filh2,B,0);
	filh2.resize(0); filh.resize(0);
	// feasible frequency range
	int numopt       = p.numopt;
	RowVector2d flim = p.flim; 
	RowVectorXi pfreq;       
	find_and(cf,flim(0),flim(1),pfreq); 
    MatrixXd  ff  = MatrixXd::Zero(nframes,numopt); 
	MatrixXd  amp = MatrixXd::Zero(nframes,numopt);
	
	VectorXi  pos;   VectorXi  pos_temp;     
	VectorXd peak;   VectorXd  peak_temp; 
	VectorXi  ind;
	RowVectorXd posff;

	RowVectorXd B_temp(pfreq.size());
	for(int ii= 1; ii<= nframes;ii++)
	{
		B_temp.resize(pfreq.size());
		for(int jj=0;jj<pfreq.size();jj++)
		    B_temp(jj) = B(ii-1,pfreq(jj)-1);
		findpeaks(B_temp,5.0/(cf(pfreq(1)-1)-cf(pfreq(0)-1)),pos_temp,peak_temp);
	    if(pos_temp.size() > 0)
		{
			peak.resize(peak_temp.size());    ind.resize(peak_temp.size());
			Sort_descend(peak_temp,peak,ind); ind = ind.array() + 1; 
			pos.resize(pos_temp.size());
			for(int nn=0;nn<pos_temp.size();nn++)
				pos(nn) = pos_temp(ind(nn)-1);
			posff.resize(pos.size());
			for(int nn=0;nn<pos.size();nn++)
				posff(nn) = cf(pfreq(pos(nn)-1)-1);
			int fin = min(numopt,posff.size());
			ff.block(ii-1,0,1,fin)  = posff.head(fin);
		    amp.block(ii-1,0,1,fin) = peak.adjoint().head(fin);
		}
	}
	B.resize(0,0);   pos.resize(0);       pos_temp.resize(0);
	peak.resize(0);  peak_temp.resize(0); ind.resize(0);
	posff.resize(0); B_temp.resize(0);
	/* Probability of the frame of being voiced */
	VectorXd pow(OO.rows()),amp_sum(OO.rows());
	for(int ii=0;ii<OO.rows();ii++){
		pow(ii) = OO.row(ii).mean() * (1E-6);
	    amp_sum(ii) = amp.row(ii).sum()/(pow(ii)*(1E9));
	}
	MatrixXd vuvfea1(pow.size(),2);
		vuvfea1 << pow,amp_sum;
    MatrixXd vuvfea2 = vuvfea1;
	OO.resize(0,0); pow.resize(0); amp_sum.resize(0);
	/* Gauss mix up */
	VectorXd pru, prv;
	gaussmixp(vuvfea1,m_u,v_u,w_u,pru,0);    //0,1-for convenient in function 
	gaussmixp(vuvfea2,m_v,v_v,w_v,prv,1);
	vuvfea1.resize(0,0); 
	vuvfea2.resize(0,0);
	//VectorXd pv = (pru-prv).array().exp()+1;
	pv = (pru-prv).array().exp()+1;
	pv = pv.array().inverse();   pru.resize(0); prv.resize(0);  
	/* Dynamic programming */
	RowVector6d w = p.w;
	MatrixXd  camp(amp.rows(),amp.cols());
	VectorXd  amp_max(amp.rows());
	for(int ii=0;ii<amp.rows();ii++)
		amp_max(ii) = amp.row(ii).maxCoeff();
	for(int ii=0;ii<amp.cols();ii++)  
		camp.col(ii) = -1*amp.col(ii).array()/amp_max.array();
	amp_max.resize(0);
	camp = (amp.array() == 0).select(w(5)*MatrixXd::Ones(amp.rows(),amp.cols()) ,camp); 
    amp.resize(0,0);
	
	double tmf  = p.tmf;
	int	   inmf	= int(tmf/txinc + 0.5);
	MatrixXd cost   = MatrixXd::Zero(nframes,numopt);
	MatrixXd prev   = cost;
	VectorXd medfx = VectorXd::Zero(nframes);
	double   dffact = 2.0/txinc;

	cost.row(0) = w(0)*camp.row(0);
	VectorXd  fpos  = ff.block(0,0,min(inmf,ff.rows()),1);
	VectorXd  med;
	get_vector(fpos,pv,inmf,0.6,med,1);
	double mf = median(med);
	if(med.size() == 0)                               //wrong in here,before solve this, imfilter is big problem
	{
		get_vector(fpos,pv,inmf,0.5,med,1);
		mf = median(med);
		if(med.size() ==0)
		{ 
			get_vector(fpos,pv,inmf,0.4,med,1); 
			mf = median(med);
			if(med.size() == 0)
			{
				get_vector(fpos,pv,inmf,0.3,med,1); 
				mf = median(med);
				if( med.size() == 0)
					mf = 0;
			}
		}
	}
	medfx(0) = mf;
	MatrixXd A1(numopt,numopt);      MatrixXd A2(numopt,numopt);
	MatrixXd df(numopt,numopt);      MatrixXd costdf(numopt,numopt);
	MatrixXd costpos(numopt,numopt); RowVector3d costf;
	for(int i= 2; i <= nframes; i++)
	{
		if(i>inmf)
		{
			fpos.resize(inmf+1);
		    fpos = ff.block(i-inmf-1,0,inmf+1,1);
			get_vector(fpos,pv,inmf,0.6,med,2); 
			mf   =  median(med);
			if(med.size() == 0)
			{
				get_vector(fpos,pv,inmf,0.5,med,2); 
				mf   =  median(med);
			    if(med.size() == 0)
				{
					get_vector(fpos,pv,inmf,0.4,med,2); 
					mf   =  median(med);
					if(med.size() == 0)
					{
						get_vector(fpos,pv,inmf,0.3,med,2); 
						mf   =  median(med);
						if(med.size() == 0)
							mf = 0;
					}
				}
			}
		}
		medfx(i-1) = mf;
		A1 << ff.row(i-1),ff.row(i-1),ff.row(i-1);
		A2 << ff.row(i-2),ff.row(i-2),ff.row(i-2);
 		df = dffact * (A1.transpose()-A2).array()/(A1.transpose()+A2).array();
		
		df = (df.array()-w(5)).square();
		costdf = df.cwiseMin(MatrixXd::Ones(numopt,numopt) * w(3))*w(2);
		
		if(med.size() == 0)
			costf << 0,0,0;
		else
			costf = (ff.row(i-1).array()-mf).abs()/mf;
		
		costpos.resize(numopt,numopt);
		costpos << cost.row(i-2),cost.row(i-2),cost.row(i-2);
		costpos = costpos + costdf; 
		MatrixXf::Index MaxIndex; 
		for(int nn=0;nn<cost.cols();nn++)
		{
			cost(i-1,nn) = costpos.row(nn).minCoeff(&MaxIndex);
		    prev(i-1,nn) =  MaxIndex+1;
		}
		cost.row(i-1) = cost.row(i-1) + w(1)*costf + w(0)*camp.row(i-1);   
	}
	// tracebank
	//VectorXd fx           = VectorXd::Zero(nframes);
	fx.resize(nframes);
			 fx           = VectorXd::Zero(nframes);
	VectorXd best		  = VectorXd::Zero(nframes);
	double cost_min		  = cost.row(cost.rows()-1).minCoeff();
	RowVectorXd cost_nose = cost.row(cost.rows()-1);
	RowVectorXi nose;
    find(cost_nose,cost_min,0,nose); 
	best(best.size()-1) = nose(0);
	fx(fx.size()-1) = ff(ff.rows()-1,best(best.size()-1)-1);
	for(int ii = nframes;ii>=2;ii--)
	{
		best(ii-2) = prev(ii-1,best(ii-1)-1);
	    fx(ii-2) = ff(ii-2,best(ii-2)-1);
	}
}
/* get the vector to calcumlte the median value */
void VAD::get_vector(VectorXd& fpos_,VectorXd &pv_,double inmf_,double flag,VectorXd &a2,int type)
{
	VectorXd  a1;
	if(type == 1){
		a1 = pv_.head(min(inmf_,pv_.size()));
	}
	if(type == 2){
		a1 = pv_.head(inmf_);
	}
	VectorXi  a1_01 = (a1.array()> flag).select(VectorXi::Ones(a1.size()),VectorXi::Zero(a1.size())); 
	a2.resize(a1_01.sum());
	int index = 0;
	for(int ii=0;ii<a1.size();ii++)
	{
		if(a1(ii) > flag)  
			a2(index++) = fpos_(ii);
	}
}
/* gaussmixp */
void VAD::gaussmixp(MatrixXd& y,MatrixXd& m,MatrixXd& v,Vector6d& w,VectorXd& lp,int flag)
{
	int		  n  = y.rows();
	int		  q  = 2;
	int		  k  = 6;
	double    fv = 1;
	int       memsize = 50000000;
			  lp = VectorXd::Zero(n);
	VectorXi  wk = VectorXi::Ones(k);

	RowVector3i  lix;   lix  << 1,2,4;
	Matrix2i     cix;   cix  << 1,2,1,2; 
	Matrix2i     rix;   rix  << 1,1,2,2;
	Matrix2i     lixi;  lixi << 1,2,2,3;
	MatrixXd     v1(6,3); 
	if(flag == 0){
       v1 << 0.01413822,   0.003357913,  0.01786169,  0.0009377269, 0.0006220489, 0.03422057,
			 0.1233703,    0.004299293,  0.007660504, 0.01779449,   0.002078821,  0.001605052,
			 1.110173,     0.00718649,   0.005734435, 0.5477135,   -0.00182316,   0.05659796;
	}
	if(flag == 1){
       v1 << 0.06181574,   0.002950501,  0.004528442, 0.2946077,    0.01433284,   0.02684239,
		     2.508473,    -0.03310555,   0.1098579,   14.17252,    -0.09009174,   0.07989255,
		     0.5834894,   -0.07854027,   0.1108958,   0.05978017,   0.005528601,  0.1309329;
	}

	int  nb  = min(n,max(1,floor(double(memsize)/(24*q*k)))); 
	int	 nl	 = ceil( double(n) / double(nb) );
	int	 jx0 = n - (nl - 1) * nb;
	RowVectorXi wnb = RowVectorXi::Ones(nb);
	RowVectorXi wnj = RowVectorXi::Ones(jx0);

	MatrixXd vi  = MatrixXd::Zero(q*k,q);
	VectorXd vim = VectorXd::Zero(q*k,1);
	VectorXd mtk = vim;
	VectorXd lvm = VectorXd::Zero(k,1);
	VectorXd wpk(12); wpk << 1,2,1,2,1,2,1,2,1,2,1,2;
	
	Matrix2d uvk,dvk_diag,veig;
	Matrix2d vik;                       //get the diagnoal value of dvk
	Vector2d dvk;					
	for(int ik = 1; ik<=k;ik++)
	{
		veig << v1(ik-1,lixi(0,0)-1),v1(ik-1,lixi(0,1)-1),v1(ik-1,lixi(1,0)-1),v1(ik-1,lixi(1,1)-1); 
        eig(veig,uvk,dvk_diag);  
		dvk(0)  = dvk_diag(0,0);  dvk(1)=dvk_diag(1,1);
		dvk_diag(0,1) = 0;        dvk_diag(0,0) = 1/dvk_diag(0,0); 
		dvk_diag(1,0) = 0;        dvk_diag(1,1) = 1/dvk_diag(1,1); 	
	    vik = -0.5*uvk * dvk_diag *uvk.adjoint();
		vi.block((ik-1)*q,0,q,q) = vik;
		vim.segment((ik-1)*q,q)  = vik * (m.row(ik-1)).adjoint();
		mtk.segment((ik-1)*q,q)  = (m.row(ik-1)).adjoint();
		lvm(ik-1) = log(w(ik-1)) - 0.5*(log(dvk(0))+log(dvk(1)));
	}

	int  jx  = jx0;
	MatrixXd xii       = y.transpose();
	MatrixXd xii_temp(vi.rows(),xii.cols()); 
	for(int ii=0;ii<xii_temp.rows();ii++)
		xii_temp.row(ii) = xii.row(wpk(ii)-1);
	MatrixXd vim_temp(vi.rows(),xii.cols()); 
	MatrixXd mtk_temp(vi.rows(),xii.cols());
	for(int ii=0; ii<vim_temp.cols();ii++){
		vim_temp.col(ii) = vim;
		mtk_temp.col(ii) = mtk;
	}
	MatrixXd A  = (vi*xii-vim_temp).array()*(xii_temp-mtk_temp).array();  
	xii.resize(0,0); xii_temp.resize(0,0); vim_temp.resize(0,0);  mtk_temp.resize(0,0);
	A.resize(q,jx*k);
	MatrixXd A1 = A.colwise().sum();  A.resize(0,0);
	A1.resize(k,jx);
	MatrixXd lvm_temp(lvm.rows(),wnj.cols());
	for(int ii=0;ii<lvm_temp.cols();ii++) 
		lvm_temp.col(ii) = lvm; 
	MatrixXd py = A1 + lvm_temp; A1.resize(0,0); lvm_temp.resize(0,0);
	RowVectorXd mx(py.cols());
	for(int ii=0; ii < py.cols();ii++) 
		mx(ii) = py.col(ii).maxCoeff();
	MatrixXd px(py.rows(),py.cols());
	for(int ii=0;ii<px.rows();ii++)  
		px.row(ii) = py.row(ii)-mx;
	px = px.array().exp();
	RowVectorXd ps = px.colwise().sum();
	lp.resize(ps.size());  
	lp = ps.array().log() + mx.array();
	// the last loop, not exist in old edition
	for(int il=2; il <= nl; il++)
	{
		cout << "Surprise! run in the gaussmixp, il loop " << endl;
		int  ix = jx+1; jx = jx+nb;
		RowVectorXi i2(jx-ix+1);
		for(int jj=0;jj<i2.size();jj++) i2(jj) = ix+jj;
		MatrixXd xii_L(y.cols(),i2.size());
		for(int jj=0;jj<i2.size();jj++)
			xii_L.col(jj) = y.row(i2(jj)-1);
			
		MatrixXd xii_temp_L(vi.rows(),xii.cols()); 
		for(int ii=0;ii<xii_temp_L.rows();ii++)
			xii_temp_L.row(ii) = xii_L.row(wpk(ii)-1);
		MatrixXd vim_temp_L(vi.rows(),xii_L.cols()); 
		MatrixXd mtk_temp_L(vi.rows(),xii_L.cols());
		for(int ii=0; ii<vim_temp_L.cols();ii++){
			vim_temp_L.col(ii) = vim;
			mtk_temp_L.col(ii) = mtk;
		}
		MatrixXd A_L  = (vi*xii_L-vim_temp_L).array()*(xii_temp_L-mtk_temp_L).array();  
		xii_L.resize(0,0);       xii_temp_L.resize(0,0); 
		vim_temp_L.resize(0,0);  mtk_temp_L.resize(0,0);
		A_L.resize(q,nb*k);
		MatrixXd A1_L = A_L.colwise().sum();  A_L.resize(0,0);
		A1_L.resize(k,nb);
		MatrixXd lvm_temp_L(lvm.rows(),wnb.cols());
		for(int ii=0;ii<lvm_temp_L.cols();ii++) 
			lvm_temp_L.col(ii) = lvm; 
		MatrixXd py_L = A1_L + lvm_temp_L; A1_L.resize(0,0); lvm_temp_L.resize(0,0);

		RowVectorXd mx_L(py_L.cols());
		for(int ii=0; ii < py_L.cols();ii++) 
			mx_L(ii) = py_L.col(ii).maxCoeff();
		MatrixXd px_L(py_L.rows(),py_L.cols());
		for(int ii=0;ii<px_L.rows();ii++)  
			px_L.row(ii) = py_L.row(ii)-mx_L;
		px_L = px_L.array().exp();
		RowVectorXd ps_L = px_L.colwise().sum();
		lp.resize(ps_L.size());  
		lp = ps_L.array().log() + mx_L.array();
	}
	lp = lp.array()-0.5*q*log(2*PI);
}

/* find peaks */
void VAD::findpeaks(RowVectorXd& x,double w,VectorXi& k,VectorXd& v)
{
   int nx = x.size();
   RowVectorXd dx(x.size()-1);
	  dx << x.tail(x.size()-1) - x.head(x.size()-1);
   RowVectorXi r;  find(dx,0,1,r);
   RowVectorXi f;  find(dx,0,-1,f);
   if(r.size()>0 && f.size()>0)
   {
	   RowVectorXi dr = r;
	   dr.tail(dr.size()-1) = r.tail(r.size()-1)-r.head(r.size()-1);
	   VectorXi rc = VectorXi::Ones(nx);
	   for(int ii=0;ii<r.size();ii++) rc(r(ii)) = 1-dr(ii); 
	   rc(0) = 0;
       VectorXi rs(rc.size());
	   int temp = 0;
	   for(int ii=0;ii<rs.size();ii++){temp=temp+rc(ii);rs(ii)=temp;}
	   
	   RowVectorXi df = f;
	   df.tail(df.size()-1) = f.tail(f.size()-1)-f.head(f.size()-1);
	   VectorXi fc = VectorXi::Ones(nx);
	   for(int ii=0;ii<f.size();ii++) fc(f(ii)) = 1-df(ii); 
	   fc(0) = 0; temp = 0;
	   VectorXi fs(fc.size());
	   for(int ii=0;ii<fs.size();ii++){temp=temp+fc(ii);fs(ii)=temp;}

	   VectorXi rp = VectorXi::Ones(nx); rp  = -1*rp;
	   rp(0) = dr(0)-1;  
	   rp(r(r.size()-1)) = nx-r(r.size()-1)-1;
	   for(int ii=0;ii<r.size()-1;ii++)
			rp(r(ii)) = dr(ii+1) -1;
	   VectorXi rq(rp.size());
	   temp = 0;
	   for(int ii=0;ii<rp.size();ii++){temp=temp+rp(ii);rq(ii)=temp;}

	   VectorXi fp = VectorXi::Ones(nx); fp  = -1*fp;
	   fp(0) = df(0)-1;  
	   fp(f(f.size()-1)) = nx-f(f.size()-1)-1;
	   for(int ii=0;ii<f.size()-1;ii++)
			fp(f(ii)) = df(ii+1) -1;
	   VectorXi fq(fp.size());
	   temp = 0;
	   for(int ii=0;ii<fp.size();ii++){temp=temp+fp(ii);fq(ii)=temp;}
	   
	   VectorXi k_temp = ((rs.array()<fs.array()) && (fq.array()<rq.array())
		    && ((fq.array()-rs.array())/2 == 0)).select(VectorXi::Ones(nx),VectorXi::Zero(nx)); 
	   
	   //VectorXi k(k_temp.sum());    
	   k.resize(k_temp.sum());     k_temp.resize(0);
	   int k_index = 0;
	   for(int ii=0;ii<nx;ii++){
			if((rs(ii)<fs(ii))&&(fq(ii)<rq(ii))&&((fq(ii)-rs(ii))/2 ==0))
				k(k_index++) = ii+1;
	   }
	   //VectorXd v(k.size());
	   v.resize(k.size());
	   for(int ii=0;ii<k.size();ii++)
			v(ii) = x(k(ii)-1);
	   RowVectorXi j;
	   RowVectorXi kk     = (k.tail(k.size()-1)-k.head(k.size()-1)).adjoint();
	   RowVectorXd kk1    = kk.cast<double>();
	   find(kk1,w,-10,j);  kk.resize(0); 

	   VectorXi k_copy;      VectorXd     v_copy; 
	   VectorXi k_copy_01;   RowVectorXi  v_01;      
	   VectorXi kk_copy;     RowVectorXd  kd_copy; 

	   int index = 0,j_num = 0;
	   while(j.any())
	   {
		   v_01.resize(j.size());
		   for(int ii=0;ii<j.size();ii++){
				if((v(j(ii)-1)) >= v(j(ii)))  v_01(ii) = 1;
				else                          v_01(ii) = 0;
		   }
		   j = j + v_01;
		   k_copy.resize(k.size());        k_copy = k;
		   v_copy.resize(v.size());        v_copy = v;
		   for(int ii=0;ii<j.size();ii++)
				k_copy(j(ii)-1) = 999999;
		   k_copy_01.resize(k_copy.size());
		   k_copy_01 = (k_copy.array()== 999999).select(VectorXi::Ones(k_copy.size()),VectorXi::Zero(k_copy.size()));
		   j_num = k_copy_01.sum();   

		   k.resize(k_copy.size()-j_num);  v.resize(v_copy.size()-j_num);
		   index = 0;
		   for(int ii=0;ii<k_copy.size();ii++){
			     if( k_copy(ii) != 999999){
					k(index)   = k_copy(ii);
					v(index++) = v_copy(ii);}
		   }
		   j.resize(0); 
		   kk_copy.resize(k.size()-1); 
		   kd_copy.resize(k.size()-1);
		   kk_copy = k.tail(k.size()-1)-k.head(k.size()-1);
		   kd_copy = kk_copy.cast<double>();
		   find(kd_copy,w,-10,j);
	   }
   }
}
/* STFT */
void VAD::spgrambw(VectorXd& s,VectorXd &t,RowVectorXd &f,MatrixXd &b)
{
	int          ns1   = s.size();
	double       bw    = 20, db = 40; 
	RowVector2d  fs(_fs,1.0/_fs);
	RowVector3d  fmax(0,5,4000);
	RowVector3d  tinc(0.01, MIN_DBL, MAX_DBL);
	double       flmin = 30;
	double       fnyq  = fs(0)/2.0;
	RowVector3d  fmaxu = fmax; 
	fmaxu(1) = fmax(1)*(fmaxu(2)-fmaxu(0))/(fmax(2)-fmax(0));   
	int          nfrq  = int((fmaxu(2)-fmaxu(0))/fmaxu(1))+1;
	RowVectorXd  fx(nfrq);
	linspace(fx,fmaxu(0),fmaxu(2));
	//RowVectorXd f;
	f = fx;
	int	   winlen  = floor(1.81 * fs(0) / bw ); 
	int	   win_n   = (winlen - (1 - winlen) ) / 2.0 + 1;
	RowVectorXd  win(win_n);
	linspace_cos(win,1- winlen,winlen);
	int	   ninc    = max(floor(tinc(0)*fs(0)+0.5),1);                                     /* window increment in samples */
	int	   fftlen  = pow(2,ceil((log((double) 4 * winlen)/log((double)2))));          /* enough oversampling to get good interpolation */
	win = win/sqrt(win.array().square().sum());
	int	 ix1  = max( floor( (tinc[1] - fs[1]) * fs[0] - (winlen - 3) / 2.0 + 0.5 ), 1 );      /* first sample required */
	int	 ix2  = min( ceil( (tinc[2] - fs[1]) * fs[0] + (winlen + 1) / 2 ), ns1 );             /* last sample required */

	VectorXd s_piece = s.segment(ix1-1,(ix2-ix1+1));
	MatrixXd sf; 
	//VectorXd t;
	enframe(s_piece,win,ninc,sf,t);
	s_piece.resize(0);  win.resize(0);
	t = fs(1) + (t.array() + ix1-2)/fs(0);
	MatrixXcd b1;
	rfft(sf,fftlen,b1);   sf.resize(0,0);
	b1 = b1.array() * b1.conjugate().array()*2.0/fs(0);
    //MatrixXd b;
	b = b1.real();        b1.resize(0,0);
	b.block(0,0,b.rows(),1)          = b.block(0,0,b.rows(),1) * 0.5; 
	b.block(0,b.cols()-1,b.rows(),1) = b.block(0,b.cols()-1,b.rows(),1)*0.5; 
	RowVectorXd  fb(fftlen/2);
	for(int ii=0;ii<fftlen/2;ii++)
		fb(ii) = ii * fs(0)/fftlen;
	double fftfs = fs(0);
	SparseMatrix<double> x;
	filtbankm_v1(nfrq,fftlen,fftfs,fx(0),fx(fx.size()-1),x);
    b = b*x.adjoint();
}

/* division frame */
void VAD::enframe(VectorXd &x,RowVectorXd &win,int inc,MatrixXd &f,VectorXd &t)
{
	int lw        = win.size();
	int nx        = x.size();
	VectorXd w    = win.transpose(); 
	int	nli	      = nx - lw + inc;                          /* ninc_ */
	int nf		  = fix(double(nli) / double(inc));       /* new number frames */
	int na        = nli - inc * nf;
	f.resize(nf,lw);
	VectorXi indf(nf);
	VectorXi inds(lw);
	for(int ii=0;ii<nf;ii++)
		indf(ii) = inc * ii; 
	for(int ii=0;ii<lw;ii++)
		inds(ii) = ii+1; 
	MatrixXi ind(nf,lw);
	double t0  = (1 + lw) / 2.0;	
	t.resize(nf);
	for(int ii=0;ii< nf;ii++)
	{
		ind.row(ii) = inds.transpose().array() + indf(ii);
		for(int jj=0;jj<lw;jj++)
		{
			f(ii,jj) = x(ind(ii,jj)-1);	
		}
		f.row(ii) = f.row(ii).array() * win.array();
		t(ii) = t0 + inc*ii;
	}
}
/* select the value of frequency input for the function of filtbankm */
void VAD::set_fin( RowVectorXd &fin1_,int nf_,double df_,RowVectorXd &fin)
{
	double* fin_  = NULL;
	int     nfin_ = 0;        
    if ( fin1_[0] == 0 )
	{
		nfin_	= nf_ + nf_ + 1;
		fin_	= new double[nfin_];
		for ( int ii = 0; ii < nf_; ii++ )
			fin_[ii] = (-1) * fin1_[nf_ - ii];
		for ( int jj = 0; jj < nf_ + 1; jj++ )
			fin_[nf_ + jj] = fin1_[jj];
	 }else if ( fin1_[0] <= df_ / 2.0 )
	{
		nfin_	= nf_ + nf_ + 1 + 1;
		fin_	= new double[nfin_];
		for ( int ii = 0; ii < nf_ + 1; ii++ )
			fin_[ii] = (-1) * fin1_[nf_ - ii];
		for ( int jj = 0; jj < nf_ + 1; jj++ )
			fin_[nf_ + jj + 1] = fin1_[jj];
	}else if ( fin1_[0] < df_ )
	{
		nfin_	= nf_ + 1 + nf_ + 1 + 2;
		fin_	= new double[nfin_];
		for ( int ii = 0; ii < nf_ + 1; ii++ )
			fin_[ii] = (-1) * fin1_[nf_ - ii];
		fin_[nf_ + 1]	= fin1_[0] - df_;
		fin_[nf_ + 2]	= df_ - fin1_[0];
		for ( int jj = 0; jj < nf_ + 1; jj++ )
			fin_[nf_ + jj + 1 + 2] = fin1_[jj];
	}else if ( fin1_[0] == df_ )
	{
		nfin_	= nf_ + nf_ + 1 + 1 + 1;
		fin_	= new double[nfin_];
		for ( int ii = 0; ii < nf_ + 1; ii++ )
			fin_[ii] = (-1) * fin1_[nf_ - ii];
		fin_[nf_ + 1] = 0;
		for ( int jj = 0; jj < nf_ + 1; jj++ )
			fin_[nf_ + jj + 1 + 1] = fin1_[jj];
	}else  {
		nfin_	= nf_ + 1 + nf_ + 1 + 2;
		fin_	= new double[nfin_];
		for ( int ii = 0; ii < nf_ + 1; ii++ )
			fin_[ii] = (-1) * fin1_[nf_ - ii];
		fin_[nf_ + 1]	= df_ - fin1_[0];
		fin_[nf_ + 2]	= fin1_[0] - df_;
		for ( int jj = 0; jj < nf_ + 1; jj++ )
			fin_[nf_ + jj + 1 + 2] = fin1_[jj];
	}
	Map<RowVectorXd> fin_temp(fin_,nfin_);                 //The same address, so can't delete fin_ 
	fin = fin_temp;    
	delete []fin_;   fin_ = NULL;
}
void VAD::stdspectrum(RowVectorXd &cf_,RowVectorXd &b)
{
	int cf_len = cf_.size();
	double	sb[7] = { 1.972049319295900E7,
					  8.911672888761445E11,
					  1.098866534797396E17,
					  2.250976456705804E21,
					  0,
					  0,
					  0 };
	double	sa[8] = { 1,
					  1.717064265336435E5,
					  8.446053794016400E9,
					  7.297481118829316E14,
					  2.255535195816080E18,
					  6.242720421841964E21,
					  3.461432392738979E24,
					  2.535193354542586E27 };
	compx	*hb	= new compx[cf_len]();
	compx	*ha	= new compx[cf_len]();
	compx	*h	= new compx[cf_len]();
	b.resize(cf_len);
	for ( int ii = 0; ii < cf_len; ii++ )
	{
		hb[ii].real = -1 * sb[0] * pow( 2 * PI * cf_[ii], 6 ) + sb[2] * pow( 2 * PI * cf_[ii], 4 )
			      - sb[4] * pow( 2 * PI * cf_[ii], 2 ) + sb[6];
		hb[ii].imag = sb[1] * pow( 2 * PI * cf_[ii], 5 ) - sb[3] * pow( 2 * PI * cf_[ii], 3 )
			      + sb[5] * 2 * PI * cf_[ii];
		ha[ii].real = -1 * sa[1] * pow( 2 * PI * cf_[ii], 6 ) + sa[3] * pow( 2 * PI * cf_[ii], 4 )
			      - sa[5] * pow( 2 * PI * cf_[ii], 2 ) + sa[7];
		ha[ii].imag = -1 * sa[0] * pow( 2 * PI * cf_[ii], 7 ) + sa[2] * pow( 2 * PI * cf_[ii], 5 )
			      - sa[4] * pow( 2 * PI * cf_[ii], 3 ) + sa[6] * (2 * PI * cf_[ii]);
		h[ii].real = (hb[ii].real * ha[ii].real + hb[ii].imag * ha[ii].imag) /
			     (ha[ii].real * ha[ii].real + ha[ii].imag * ha[ii].imag);
		h[ii].imag = (hb[ii].imag * ha[ii].real - hb[ii].real * ha[ii].imag) /
			     (ha[ii].real * ha[ii].real + ha[ii].imag * ha[ii].imag);
		b[ii] = h[ii].real * h[ii].real + h[ii].imag * h[ii].imag;
	}
	delete[]hb;  hb = NULL; delete[]ha;  ha = NULL; delete[]h;   h  = NULL;
}
void VAD::smooth(RowVectorXd &x,int n) /* eltass = smooth(eltass,p.mpsmooth); */
{
	int          nx = x.size();
	RowVectorXd  c(nx);
	double temp = 0;
	for(int ii=0;ii<nx;ii++){
		temp	= temp + x(ii); c(ii)	= temp;}
	n = 201;                           /* n = p.mpsmooth */
	int N[101] = {};
	for (int i = 0; i < 101; i++ )
	{
		N[i] = 2 * i + 1;
		x[i] = c[N[i] - 1] / N[i];
	}
	int n1 = nx - (n + 1) + 1;
	for (int i = 101; i < 101 + n1; i++ )
	{
		x[i] = (c[(n + 1 - 1) + (i - 101)] - c[i - 101]) / double(n);
	}
	for (int i = 101 + n1; i < nx; i++ )
	{
		x[i] = (c[nx - 1] - c[(nx - n + 2 - 1) + 2 * (i - 101 - n1)]) / double(n - 2 - 2 * (i - 101 - n1) );
	}
}