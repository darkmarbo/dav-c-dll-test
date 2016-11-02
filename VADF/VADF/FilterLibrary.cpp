/*-------------------------------------------------------------------*/
/**  Name: fufill all filter banks	                                **/
/*-------------------------------------------------------------------*/
#include "VAD.h"

void VAD::imfilter(MatrixXd &OO_,VectorXd& B,MatrixXd& AF_,int flag)
{
	/*data copy */
	double *OO1  = new double[OO_.rows() * OO_.cols()]();
	double *AF1  = new double[OO_.rows() * OO_.cols()]();
	double * B1  = new double[B.size()](); 
	for(int ii=0;ii<OO_.rows();ii++)
	{
		for(int jj=0;jj<OO_.cols();jj++)
		{
			OO1[ii*OO_.cols() + jj] = OO_(ii,jj);
		}
	}
	for(int ii= 0 ;ii<B.size();ii++)
		B1[ii] = B(ii); 
	int OO_row = OO_.rows();
	int OO_col = OO_.cols();
	int B_len  = B.size();
	/* copy over */

	int	R_core	= 0;
	int	L_core	= 0;
	int	core	= 0;
	if ( B_len % 2 == 0 )
		core = B_len / 2;
	else
		core = (B_len + 1) / 2;
	L_core	= core - 1;
	R_core	= B_len - core;
	int	OO_L	= 0;
	int	OO_R	= 0;
	int	L	= 0;    /* min(L_core,OO_L) */
	int	R	= 0;    /* min(R_core,OO_R) */
	if ( flag == 0 )
	{
		for ( int ii = 0; ii < OO_row; ii++ )
		{
			for ( int jj = 0; jj < OO_col; jj++ )
			{
				OO_L	= (jj + 1) - 1;
				OO_R	= OO_col - (jj + 1);
				L	= min( L_core, OO_L );
				R	= min( R_core, OO_R );
				for ( int nn = 0; nn < R + L + 1; nn++ )
				{
					AF1[ii * OO_col + jj] = AF1[ii * OO_col + jj] + OO1[ii * OO_col + (jj - L) + nn]
							       * B1[(core - L - 1) + nn];
				}
			}
		}
	}

	if ( flag == 1 )
	{
		for ( int ii = 0; ii < OO_row; ii++ )
		{
			for ( int jj = 0; jj < OO_col; jj++ )
			{
				OO_L	= (ii + 1) - 1;
				OO_R	= OO_row - (ii + 1);
				L	= min( L_core, OO_L );
				R	= min( R_core, OO_R );
				for ( int nn = 0; nn < R + L + 1; nn++ )
				{
					AF1[ii * OO_col + jj] = AF1[ii * OO_col + jj] + OO1[(ii - L + nn) * OO_col + jj]
							       * B1[(core - L - 1) + nn];
				}
			}
		}
	}
	for(int ii=0;ii<OO_.rows();ii++)
	{
		for(int jj=0;jj<OO_.cols();jj++)
		{
			AF_(ii,jj) = AF1[ii*OO_.cols() + jj];
		}
	}
	delete []OO1;   OO1 = NULL;
	delete []AF1;   AF1 = NULL;
	delete []B1;    B1  = NULL;
}  
// design filtbankm_v2
void VAD::filtbankm_v2(int p,int n,double fs,double fl,double fh,SparseMatrix<double> &x,RowVectorXd &cf) 
{
	double	f11 = 0;     //not the same as fl param
	int     nf  = 1+floor(double(n)/2.0);
	double	df  = fs/n;
	RowVectorXd fin0(nf);
	for(int ii=0;ii<nf;ii++)
		fin0(ii) = f11 + ii*df;
	RowVector2d  mflh(fl,fh);
	mflh(0)=log10(mflh(0)); 
	mflh(1)=log10(mflh(1));
	double melrng = mflh(0)*(-1)+mflh(1);
	double melinc = melrng/(p+1);
    //RowVectorXd  cf(p+2);
	cf.resize(p+2);
	for(int ii=0;ii<p+2;ii++){
		cf(ii) = mflh(0)+ii*melinc;
	    if(ii>0) cf(ii) = max(cf(ii),0); 
	}
	RowVectorXd mb(cf.size());
	for(int ii=0;ii<cf.size();ii++) 
		mb(ii) = pow(10,cf(ii));

	RowVectorXd fin1(nf+1); 
	fin1.head(nf) = fin0;
	fin1(nf)      = fin1(nf-1) + df;
	RowVectorXd fin; 
	set_fin(fin1,nf,df,fin);  fin1.resize(0);
	int  nfin  = fin.size();

	RowVectorXd  fout = mb;
	int  mfout = fout.size();
	RowVectorXd  gout = fout.segment(2,mfout-2) - fout.segment(0,mfout-2); 
	RowVectorXd  vv   = (gout.array() == 0).select(RowVectorXd::Ones(gout.size()),RowVectorXd::Zero(gout.size()));
	gout = 2*(gout+vv).array().inverse();  vv.resize(0); 

    RowVectorXd  gin  = RowVectorXd::Ones(nfin-2);

	RowVectorXd foutin(fout.size()+fin.size());
	     foutin << fout,fin;
    int  nfall = foutin.size();
	RowVectorXd wleft(1+(mfout-1)+1+(nfin-1));
	     wleft  << 0,fout.segment(1,mfout-1)-fout.segment(0,mfout-1),0,fin.segment(1,nfin-1) - fin.segment(0,nfin-1);
    RowVectorXd wright(wleft.size());
	     wright << wleft.segment(1,wleft.size()-1),0;
	RowVectorXd ffact(3+gout.size()+ min(nf,nfin-nf-2)+max(nfin-2*nf-2,0)+ (nfin-2)-(nfin-nf-1)+2);
	     ffact  << 0,gout,0,0,gin.head(min(nf,nfin-nf-2)),RowVectorXd::Zero(max(nfin-2*nf-2,0)),gin.segment(nfin-nf-2,(nfin-2)-(nfin-nf-1)+1),0; 
    
	RowVectorXi ifall(foutin.size());                //be careful  with this vector
	RowVectorXd fall(foutin.size());
	sort_up(foutin,fall,ifall);                      //排序后数组位置比matlab小1，也即从0开始 
	ifall.array() = ifall.array()+1;                          // now,the same with matlab
	RowVectorXi jfall = RowVectorXi::Zero(nfall);
	RowVectorXi infall(nfall);
	for(int ii=0;ii<nfall;ii++){
		infall(ii) = ii+1; 
		jfall(ifall(ii)-1) = infall(ii);
	}
	int njj1 = max(jfall(0),jfall(mfout))-2;
	int njj2 = min(jfall(mfout-1),jfall(nfall-1))+2;
	RowVectorXi ifall_pic(njj1+(nfall-njj2+1));
	    ifall_pic << ifall.head(njj1),ifall.segment(njj2-1,nfall-njj2+1);
	for(int ii=0;ii<njj1+(nfall-njj2+1);ii++)
		ffact(ifall_pic(ii)-1) = 0;
	ifall_pic.resize(0);
	
	RowVectorXi nxto = (ifall.array() <= mfout).select(RowVectorXi::Ones(ifall.size()),RowVectorXi::Zero(ifall.size()));
	RowVectorXi nxti = (ifall.array() >  mfout).select(RowVectorXi::Ones(ifall.size()),RowVectorXi::Zero(ifall.size()));
	RowVectorXi nxtr(nxti.size()); 
	int	temp1 = 0,temp2	= 0;	
	for (int ii=0;ii<ifall.size();ii++ ){
		temp1    = temp1 + nxto[ii];  nxto[ii]	= temp1;
		temp2    = temp2 + nxti[ii];  nxti[ii]	= temp2;
		nxtr(ii) = min(nxti(ii)+1+mfout,nfall);
	}
	ifall_pic = (ifall.array() > mfout).select(RowVectorXi::Ones(ifall.size()),RowVectorXi::Zero(ifall.size()));
	for(int ii=0;ii<ifall.size();ii++){
		if(ifall_pic(ii) == 1)
			nxtr(ii) = 1 + nxto(ii);
	} 
	ifall_pic.resize(0);
	RowVectorXi nxtr_bk = nxtr;
	for(int ii=0;ii<nxtr.size();ii++)
		nxtr(ii) = nxtr_bk(jfall(ii)-1);
	nxtr_bk.resize(0);
	/* part 1  */
	RowVectorXi msk0 =(ffact.array()>0).select(RowVectorXi::Ones(ffact.size()),RowVectorXi::Zero(ffact.size()));
	//msk.resize(msk0.size());
	RowVectorXi msk(msk0.size());
	int flag = 0;
	for(int ii=0;ii<msk0.size();ii++){
		if( ffact(nxtr(ii)-1) > 0) flag = 1;
		else                       flag = 0;
		msk(ii) = msk0(ii) & flag;
	}
	RowVectorXi ix1(msk.sum());
	RowVectorXi jx1(msk.sum());
	int index = 0;
	for(int ii=0;ii<msk.size();ii++){
		if(msk(ii) == 1){
			ix1(index)   = infall(ii);
			jx1(index++) = nxtr(ii);}
	}
	RowVectorXd vfgx(ix1.size());
	RowVectorXd yx(vfgx.size());
	for(int ii=0;ii<ix1.size();ii++){
		vfgx(ii) = foutin(ix1(ii)-1) - foutin(jx1(ii)-2);
        yx(ii)   = min(wleft(ix1(ii)-1),vfgx(ii));
	}
	RowVectorXd yx_01 = (yx.array() == 0).select(RowVectorXd::Ones(yx.size()),RowVectorXd::Zero(yx.size()));
	RowVectorXd wx1(vfgx.size());
	for(int ii=0;ii<vfgx.size();ii++){
		wx1(ii) = ffact(ix1(ii)-1)*ffact(jx1(ii)-1)*yx(ii)*(wleft(ix1(ii)-1)*vfgx(ii)-yx(ii)*(0.5*(wleft(ix1(ii)-1)
			      + vfgx(ii))-yx(ii)/3.0))/(wleft(ix1(ii)-1)*wleft(jx1(ii)-1) + yx_01(ii)); 
	}
	yx_01.resize(0);
	/* part2 */
	RowVectorXi nxtu(nxtr.size());
	    nxtu << nxtr.segment(1,nxtr.size()-1).array()-1,0;
	for(int ii=0;ii<nxtr.size();ii++) nxtu(ii) = max(nxtu(ii),1);
	flag = 0;
	for(int ii=0;ii<msk0.size();ii++){
		if( ffact(nxtu(ii)-1) > 0) flag = 1;
		else                       flag = 0;
		msk(ii) = msk0(ii) & flag;
	}
	RowVectorXi ix2(msk.sum());
	RowVectorXi jx2(msk.sum());
	index = 0;
	for(int ii=0;ii<msk.size();ii++){
		if(msk(ii) == 1){
			ix2(index)   = infall(ii);
			jx2(index++) = nxtu(ii);}
	}
	vfgx.resize(ix2.size());
	yx.resize(ix2.size());
	for(int ii=0;ii<ix2.size();ii++){
		vfgx(ii) = foutin(ix2(ii)) - foutin(jx2(ii)-1);
        yx(ii)   = min(wright(ix2(ii)-1),vfgx(ii));
		if(foutin(jx2(ii)) < foutin(ix2(ii)))
			yx(ii) = 0;
	}
	yx_01.resize(yx.size());
	yx_01 = (yx.array() == 0).select(RowVectorXd::Ones(yx.size()),RowVectorXd::Zero(yx.size()));
    RowVectorXd wx2(vfgx.size());
	for(int ii=0;ii<vfgx.size(); ii++ ){
		wx2(ii) = ffact(ix2(ii)- 1) * ffact(jx2(ii)- 1) * yx(ii) * yx(ii)*(0.5 * (wright(jx2(ii)-1)
			     - vfgx(ii)) + yx(ii)/3.0)/(wright(ix2(ii)-1) * wright(jx2(ii)-1)+yx_01(ii));
	}
	/* part3 */
	nxtu.resize(nxtr.size());
	for(int ii=0;ii<nxtr.size();ii++) nxtu(ii) = max(nxtr(ii)-1,1);
	flag = 0;
	for(int ii=0;ii<msk0.size();ii++){
		if( ffact(nxtu(ii)-1) > 0) flag = 1;
		else                       flag = 0;
		msk(ii) = msk0(ii) & flag;
	}
	RowVectorXi ix3(msk.sum());
	RowVectorXi jx3(msk.sum());
	index = 0;
	for(int ii=0;ii<msk.size();ii++){
		if(msk(ii) == 1){
			ix3(index)   = infall(ii);
			jx3(index++) = nxtu(ii);}
	}
	vfgx.resize(ix3.size());
	yx.resize(ix3.size());
	for(int ii=0;ii<ix3.size();ii++){
		vfgx(ii) = foutin(ix3(ii)-1) - foutin(jx3(ii)-1);
        yx(ii)   = min(wleft(ix3(ii)-1),vfgx(ii));
		if(foutin(jx3(ii)) < foutin(ix3(ii)-1))
			yx(ii) = 0;
	}
	yx_01.resize(yx.size());
	yx_01 = (yx.array() == 0).select(RowVectorXd::Ones(yx.size()),RowVectorXd::Zero(yx.size()));
	RowVectorXd wx3(vfgx.size());
	for(int ii=0;ii<vfgx.size(); ii++ ){
		wx3(ii) = ffact(ix3(ii)-1) * ffact(jx3(ii)-1)*yx(ii) * (wleft(ix3(ii)-1) * (wright(jx3(ii)-1)- vfgx(ii)) + yx(ii)*
			   (0.5 * (wleft(ix3(ii)-1) - wright(jx3(ii)-1) + vfgx(ii))-yx(ii)/3.0))/(wleft(ix3(ii)-1) * wright(jx3(ii)-1)+yx_01(ii)); 
	}
	/*  part4 */
	nxtu.resize(nxtr.size());
	    nxtu << nxtr.segment(1,nxtr.size()-1),1;
	flag = 0;
	for(int ii=0;ii<msk0.size();ii++){
		if( ffact(nxtu(ii)-1) > 0) flag = 1;
		else                       flag = 0;
		msk(ii) = msk0(ii) & flag;
	}
	RowVectorXi ix4(msk.sum());
	RowVectorXi jx4(msk.sum());
	index = 0;
	for(int ii=0;ii<msk.size();ii++){
		if(msk(ii) == 1){
			ix4(index)   = infall(ii);
			jx4(index++) = nxtu(ii);}
	}
	vfgx.resize(ix4.size());
	yx.resize(ix4.size());
	for(int ii=0;ii<ix4.size();ii++){
		vfgx(ii) = foutin(ix4(ii)) - foutin(jx4(ii)-2);
        yx(ii)   = min(wright(ix4(ii)-1),vfgx(ii));
	}
	yx_01.resize(yx.size());
	yx_01 = (yx.array() == 0).select(RowVectorXd::Ones(yx.size()),RowVectorXd::Zero(yx.size()));
	RowVectorXd wx4(vfgx.size());
	for(int ii=0;ii<vfgx.size(); ii++){
		wx4[ii] = ffact(ix4(ii)-1)*ffact(jx4(ii)-1)*yx(ii)*yx(ii)*(0.5*vfgx(ii)
				  -yx(ii)/3.0)/(wright(ix4(ii)-1)*wleft(jx4(ii)-1)+yx_01(ii));
	}
	/* in sum  */
	MatrixXi iox(2,ix1.size()+ix2.size()+ix3.size()+ix4.size());
	   iox << ix1,ix2,ix3,ix4,jx1,jx2,jx3,jx4 ;
	temp1 = 0; temp2 = 0;
    for(int ii=0;ii<ix1.size()+ix2.size()+ix3.size()+ix4.size();ii++)
	{
		temp1 = iox(0,ii);  temp2 = iox(1,ii); 
		iox(0,ii) = min(temp1,temp2);
		iox(1,ii) = max(temp1,temp2);
	}
    msk.resize(iox.cols());
	msk = ((nfall+mfout)/2.0 >= iox.block(1,0,1,iox.cols()).array()).select(RowVectorXi::Ones(msk.size()),RowVectorXi::Zero(msk.size())); 
	for(int ii=0;ii<msk.size();ii++)
	{
		if(msk(ii) == 1) iox(1,ii) = (nfall+mfout+1)-iox(1,ii);
	}
	RowVectorXd S(wx1.size()+wx2.size()+wx3.size()+wx4.size());
		S << wx1,wx2,wx3,wx4;
	RowVectorXi ps1(iox.cols());
		ps1 << iox.block(0,0,1,iox.cols()).array() -1;
	RowVectorXi ps2(iox.cols());
	    ps2 << iox.block(1,0,1,iox.cols()).array()-nfall+nf+1;
	for(int ii=0;ii<ps2.size();ii++) ps2(ii) = max(ps2(ii),1);
	/* Generate SparseMatrix */
	//SparseMatrix<double> x(p,nf);
	x.resize(p,nf);
	vector<Triplet<double>> triplets;    /* ps1,ps2行列 */
	for(int ii=0;ii<S.size();ii++)
	   triplets.push_back(Triplet<double>(ps1(ii)-1,ps2(ii)-1,S(ii)));		
	x.setFromTriplets(triplets.begin(),triplets.end()) ;
	cf.resize(p);
	cf = mb.segment(1,p);         //for the filtbank 2
	VectorXd   sx(x.rows());
	for(int ii=0;ii<x.rows();ii++)
	   sx(ii) = x.block(ii,0,1,x.cols()).sum();
	VectorXd   sx_temp = (sx.array()==0).select(VectorXd::Ones(sx.size()),VectorXd::Zero(sx.size()));  
	sx = sx + sx_temp;   
	for(int jj=0;jj<x.cols();jj++)   // block 不能为左值
		x.col(jj) = x.block(0,jj,x.rows(),1).cwiseProduct(sx.cwiseInverse());
}
// design filtbankm_v1
void VAD::filtbankm_v1(int p,int n,double fs,double fl,double fh,SparseMatrix<double> &x) 
{
	double f11 = 0;
	int     nf  = 1+floor(double(n)/2.0);
	double	df  = fs/n;
	RowVectorXd fin0(nf);
	for(int ii=0;ii<nf;ii++)
		fin0(ii) = f11 + ii*df;
	RowVector2d  mflh(fl,fh);
	double melrng = mflh(0)*(-1)+mflh(1);
	double melinc = melrng/(p-1);
	mflh(0) = mflh(0) - melinc;
	mflh(1) = mflh(1) + melinc;
    RowVectorXd  cf(p+2);
	for(int ii=0;ii<p+2;ii++){
		cf(ii) = mflh(0)+ii*melinc;
	    if(ii>0) cf(ii) = max(cf(ii),0); 
	}
	RowVectorXd mb  = cf;
	RowVectorXd fin1(nf+1); 
	fin1.head(nf) = fin0;
	fin1(nf)      = fin1(nf-1) + df;
	RowVectorXd fin; 
	set_fin(fin1,nf,df,fin);  fin1.resize(0);
	int  nfin  = fin.size();

	RowVectorXd  fout = mb;
	int  mfout = fout.size();
	RowVectorXd  gout = fout.segment(2,mfout-2) - fout.segment(0,mfout-2); 
	RowVectorXd  vv   = (gout.array() == 0).select(RowVectorXd::Ones(gout.size()),RowVectorXd::Zero(gout.size()));
	gout = 2*(gout+vv).array().inverse();  vv.resize(0); 

    RowVectorXd  gin  = RowVectorXd::Ones(nfin-2);
	RowVectorXi  msk  = (fin.segment(1,nfin-2).array() == 0).select(RowVectorXi::Ones(nfin-2),RowVectorXi::Zero(nfin-2));
	for(int ii=0;ii<msk.size();ii++){
		if(msk(ii) == 0)
			gin(ii) = 0.5 * gin(ii);
	}

	RowVectorXd foutin(fout.size()+fin.size());
	     foutin << fout,fin;
    int  nfall = foutin.size();
	RowVectorXd wleft(1+(mfout-1)+1+(nfin-1));
	     wleft  << 0,fout.segment(1,mfout-1)-fout.segment(0,mfout-1),0,fin.segment(1,nfin-1) - fin.segment(0,nfin-1);
    RowVectorXd wright(wleft.size());
	     wright << wleft.segment(1,wleft.size()-1),0;
	RowVectorXd ffact(3+gout.size()+ min(nf,nfin-nf-2)+max(nfin-2*nf-2,0)+ (nfin-2)-(nfin-nf-1)+2);
	     ffact  << 0,gout,0,0,gin.head(min(nf,nfin-nf-2)),RowVectorXd::Zero(max(nfin-2*nf-2,0)),gin.segment(nfin-nf-2,(nfin-2)-(nfin-nf-1)+1),0; 
    
	RowVectorXi ifall(foutin.size());                //be careful  with this vector
	RowVectorXd fall(foutin.size());
	sort_up(foutin,fall,ifall);                      //排序后数组位置比matlab小1，也即从0开始 
	ifall.array() = ifall.array()+1;                          // now,the same with matlab
	RowVectorXi jfall = RowVectorXi::Zero(nfall);
	RowVectorXi infall(nfall);
	for(int ii=0;ii<nfall;ii++){
		infall(ii) = ii+1; 
		jfall(ifall(ii)-1) = infall(ii);
	}
	int njj1 = max(jfall(0),jfall(mfout))-2;
	int njj2 = min(jfall(mfout-1),jfall(nfall-1))+2;
	RowVectorXi ifall_pic(njj1+(nfall-njj2+1));
	    ifall_pic << ifall.head(njj1),ifall.segment(njj2-1,nfall-njj2+1);
	for(int ii=0;ii<njj1+(nfall-njj2+1);ii++)
		ffact(ifall_pic(ii)-1) = 0;
	ifall_pic.resize(0);
	
	RowVectorXi nxto = (ifall.array() <= mfout).select(RowVectorXi::Ones(ifall.size()),RowVectorXi::Zero(ifall.size()));
	RowVectorXi nxti = (ifall.array() >  mfout).select(RowVectorXi::Ones(ifall.size()),RowVectorXi::Zero(ifall.size()));
	RowVectorXi nxtr(nxti.size()); 
	int	temp1 = 0,temp2	= 0;	
	for (int ii=0;ii<ifall.size();ii++ ){
		temp1    = temp1 + nxto[ii];  nxto[ii]	= temp1;
		temp2    = temp2 + nxti[ii];  nxti[ii]	= temp2;
		nxtr(ii) = min(nxti(ii)+1+mfout,nfall);
	}
	ifall_pic = (ifall.array() > mfout).select(RowVectorXi::Ones(ifall.size()),RowVectorXi::Zero(ifall.size()));
	for(int ii=0;ii<ifall.size();ii++){
		if(ifall_pic(ii) == 1)
			nxtr(ii) = 1 + nxto(ii);
	} 
	ifall_pic.resize(0);
	RowVectorXi nxtr_bk = nxtr;
	for(int ii=0;ii<nxtr.size();ii++)
		nxtr(ii) = nxtr_bk(jfall(ii)-1);
	nxtr_bk.resize(0);
	/* part 1  */
	RowVectorXi msk0 =(ffact.array()>0).select(RowVectorXi::Ones(ffact.size()),RowVectorXi::Zero(ffact.size()));
	msk.resize(msk0.size());
	int flag = 0;
	for(int ii=0;ii<msk0.size();ii++){
		if( ffact(nxtr(ii)-1) > 0) flag = 1;
		else                       flag = 0;
		msk(ii) = msk0(ii) & flag;
	}
	RowVectorXi ix1(msk.sum());
	RowVectorXi jx1(msk.sum());
	int index = 0;
	for(int ii=0;ii<msk.size();ii++){
		if(msk(ii) == 1){
			ix1(index)   = infall(ii);
			jx1(index++) = nxtr(ii);}
	}
	RowVectorXd vfgx(ix1.size());
	RowVectorXd yx(vfgx.size());
	for(int ii=0;ii<ix1.size();ii++){
		vfgx(ii) = foutin(ix1(ii)-1) - foutin(jx1(ii)-2);
        yx(ii)   = min(wleft(ix1(ii)-1),vfgx(ii));
	}
	RowVectorXd yx_01 = (yx.array() == 0).select(RowVectorXd::Ones(yx.size()),RowVectorXd::Zero(yx.size()));
	RowVectorXd wx1(vfgx.size());
	for(int ii=0;ii<vfgx.size();ii++){
		wx1(ii) = ffact(ix1(ii)-1)*ffact(jx1(ii)-1)*yx(ii)*(wleft(ix1(ii)-1)*vfgx(ii)-yx(ii)*(0.5*(wleft(ix1(ii)-1)
			      + vfgx(ii))-yx(ii)/3.0))/(wleft(ix1(ii)-1)*wleft(jx1(ii)-1) + yx_01(ii)); 
	}
	yx_01.resize(0);
	/* part2 */
	RowVectorXi nxtu(nxtr.size());
	    nxtu << nxtr.segment(1,nxtr.size()-1).array()-1,0;
	for(int ii=0;ii<nxtr.size();ii++) nxtu(ii) = max(nxtu(ii),1);
	flag = 0;
	for(int ii=0;ii<msk0.size();ii++){
		if( ffact(nxtu(ii)-1) > 0) flag = 1;
		else                       flag = 0;
		msk(ii) = msk0(ii) & flag;
	}
	RowVectorXi ix2(msk.sum());
	RowVectorXi jx2(msk.sum());
	index = 0;
	for(int ii=0;ii<msk.size();ii++){
		if(msk(ii) == 1){
			ix2(index)   = infall(ii);
			jx2(index++) = nxtu(ii);}
	}
	vfgx.resize(ix2.size());
	yx.resize(ix2.size());
	for(int ii=0;ii<ix2.size();ii++){
		vfgx(ii) = foutin(ix2(ii)) - foutin(jx2(ii)-1);
        yx(ii)   = min(wright(ix2(ii)-1),vfgx(ii));
		if(foutin(jx2(ii)) < foutin(ix2(ii)))
			yx(ii) = 0;
	}
	yx_01.resize(yx.size());
	yx_01 = (yx.array() == 0).select(RowVectorXd::Ones(yx.size()),RowVectorXd::Zero(yx.size()));
    RowVectorXd wx2(vfgx.size());
	for(int ii=0;ii<vfgx.size(); ii++ )
	{
		wx2(ii) = ffact(ix2(ii)- 1) * ffact(jx2(ii)- 1) * yx(ii) * yx(ii)*(0.5 * (wright(jx2(ii)-1)
			     - vfgx(ii)) + yx(ii)/3.0)/(wright(ix2(ii)-1) * wright(jx2(ii)-1)+yx_01(ii));
	}
	/* part3 */
	nxtu.resize(nxtr.size());
	for(int ii=0;ii<nxtr.size();ii++) nxtu(ii) = max(nxtr(ii)-1,1);
	flag = 0;
	for(int ii=0;ii<msk0.size();ii++){
		if( ffact(nxtu(ii)-1) > 0) flag = 1;
		else                       flag = 0;
		msk(ii) = msk0(ii) & flag;
	}
	RowVectorXi ix3(msk.sum());
	RowVectorXi jx3(msk.sum());
	index = 0;
	for(int ii=0;ii<msk.size();ii++){
		if(msk(ii) == 1){
			ix3(index)   = infall(ii);
			jx3(index++) = nxtu(ii);}
	}
	vfgx.resize(ix3.size());
	yx.resize(ix3.size());
	for(int ii=0;ii<ix3.size();ii++){
		vfgx(ii) = foutin(ix3(ii)-1) - foutin(jx3(ii)-1);
        yx(ii)   = min(wleft(ix3(ii)-1),vfgx(ii));
		if(foutin(jx3(ii)) < foutin(ix3(ii)-1))
			yx(ii) = 0;
	}
	yx_01.resize(yx.size());
	yx_01 = (yx.array() == 0).select(RowVectorXd::Ones(yx.size()),RowVectorXd::Zero(yx.size()));
	RowVectorXd wx3(vfgx.size());
	for(int ii=0;ii<vfgx.size(); ii++ )
	{
		wx3(ii) = ffact(ix3(ii)-1) * ffact(jx3(ii)-1)*yx(ii) * (wleft(ix3(ii)-1) * (wright(jx3(ii)-1)- vfgx(ii)) + yx(ii)*
			   (0.5 * (wleft(ix3(ii)-1) - wright(jx3(ii)-1) + vfgx(ii))-yx(ii)/3.0))/(wleft(ix3(ii)-1) * wright(jx3(ii)-1)+yx_01(ii)); 
	}
	/*  part4 */
	nxtu.resize(nxtr.size());
	    nxtu << nxtr.segment(1,nxtr.size()-1),1;
	flag = 0;
	for(int ii=0;ii<msk0.size();ii++){
		if( ffact(nxtu(ii)-1) > 0) flag = 1;
		else                       flag = 0;
		msk(ii) = msk0(ii) & flag;
	}
	RowVectorXi ix4(msk.sum());
	RowVectorXi jx4(msk.sum());
	index = 0;
	for(int ii=0;ii<msk.size();ii++){
		if(msk(ii) == 1){
			ix4(index)   = infall(ii);
			jx4(index++) = nxtu(ii);}
	}
	vfgx.resize(ix4.size());
	yx.resize(ix4.size());
	for(int ii=0;ii<ix4.size();ii++){
		vfgx(ii) = foutin(ix4(ii)) - foutin(jx4(ii)-2);
        yx(ii)   = min(wright(ix4(ii)-1),vfgx(ii));          //wrong in here wright
	}
	yx_01.resize(yx.size());
	yx_01 = (yx.array() == 0).select(RowVectorXd::Ones(yx.size()),RowVectorXd::Zero(yx.size()));
	RowVectorXd wx4(vfgx.size());
	for(int ii=0;ii<vfgx.size(); ii++)
	{
		wx4[ii] = ffact(ix4(ii)-1)*ffact(jx4(ii)-1)*yx(ii)*yx(ii)*(0.5*vfgx(ii)
				  -yx(ii)/3.0)/(wright(ix4(ii)-1)*wleft(jx4(ii)-1)+yx_01(ii));
	}
	/* in sum  */
	MatrixXi iox(2,ix1.size()+ix2.size()+ix3.size()+ix4.size());
	   iox << ix1,ix2,ix3,ix4,jx1,jx2,jx3,jx4 ;
	temp1 = 0; temp2 = 0;
    for(int ii=0;ii<ix1.size()+ix2.size()+ix3.size()+ix4.size();ii++)
	{
		temp1 = iox(0,ii);  temp2 = iox(1,ii); 
		iox(0,ii) = min(temp1,temp2);
		iox(1,ii) = max(temp1,temp2);
	}
    msk.resize(iox.cols());
	msk = ((nfall+mfout)/2.0 >= iox.block(1,0,1,iox.cols()).array()).select(RowVectorXi::Ones(msk.size()),RowVectorXi::Zero(msk.size())); 
	for(int ii=0;ii<msk.size();ii++)
	{
		if(msk(ii) == 1) iox(1,ii) = (nfall+mfout+1)-iox(1,ii);
	}
	RowVectorXd S(wx1.size()+wx2.size()+wx3.size()+wx4.size());
		S << wx1,wx2,wx3,wx4;
	RowVectorXi ps1(iox.cols());
		ps1 << iox.block(0,0,1,iox.cols()).array() -1;
	RowVectorXi ps2(iox.cols());
	    ps2 << iox.block(1,0,1,iox.cols()).array()-nfall+nf+1;
	for(int ii=0;ii<ps2.size();ii++) ps2(ii) = max(ps2(ii),1);
	/* Generate SparseMatrix */
	//SparseMatrix<double> x(p,nf);
	x.resize(p,nf);
	vector<Triplet<double>> triplets;    /* ps1,ps2行列 */
	for(int ii=0;ii<S.size();ii++)
	   triplets.push_back(Triplet<double>(ps1(ii)-1,ps2(ii)-1,S(ii)));		
	x.setFromTriplets(triplets.begin(),triplets.end()) ;

	cf.resize(p);
	cf = mb.segment(1,p);         //for the filtbank 2
	VectorXd   sx(x.rows());
	for(int ii=0;ii<x.rows();ii++)
	   sx(ii) = x.block(ii,0,1,x.cols()).sum();
	VectorXd   sx_temp = (sx.array()==0).select(VectorXd::Ones(sx.size()),VectorXd::Zero(sx.size()));  
	sx = sx + sx_temp;   
	for(int jj=0;jj<x.cols();jj++)   // block 不能为左值
		x.col(jj) = x.block(0,jj,x.rows(),1).cwiseProduct(sx.cwiseInverse());
}