/*-------------------------------------------------------------------*/
/**  Name: fufill the core part of VAD ---VFR	                    **/
/*-------------------------------------------------------------------*/
#include "VAD.h"

void VAD::vad_snre_pv(VectorXd& dfdata,VectorXd& pv01,VectorXd& pvblk,double vadThres,VectorXd& pv_vad)
{
	int	Dexpl	= 18;
	int	Dexpr	= 18;
	VectorXd Dsmth = VectorXd::Zero(_nfr10);
	VectorXd e     = VectorXd::Zero(_nfr10);
	for(int ii=1;ii <= _nfr10;ii++)
	{
		for(int jj=1;jj <= _flen; jj++)
		{
			e[ii-1] = e[ii-1]+dfdata[(ii-1)*_fsh10+jj-1]*dfdata[(ii-1)*_fsh10+jj-1];
		}
		if (e[ii-1] <= ENERGYFLOOR)
			e[ii-1] = ENERGYFLOOR;
	}
	double	 segsnrsmth	 = 1;
	double	 sign_segsnr = 0;
	VectorXd segsnr   = VectorXd::Zero(_nfr10);
	VectorXd D        = VectorXd::Zero(_nfr10);
	VectorXd postsnr  = VectorXd::Zero(_nfr10);
	VectorXd snre_vad = VectorXd::Zero(_nfr10);
	int	     sign_pv  = 0;
	int	     nstart   = 0;
	int	     nstop    = 0;
	// all the temporary matrix in the next loop
	VectorXd datai;  VectorXd Dexp; VectorXd Dexp_temp(Dexpl+Dexpr);
	VectorXd e_temp; VectorXd eY;
	double emin;
	double snr_adj = 4;

	for(int ii=1; ii<= _nfr10;ii++)
	{
		if ((pvblk[ii-1]== 1) && (sign_pv == 0))
		{
			nstart	= ii;
			sign_pv = 1;
		}
		else if((pvblk[ii-1] == 0 || ii == _nfr10) && (sign_pv == 1)) /* fuck1 */
		{
			nstop = ii-1;
			if(ii == _nfr10 ) nstop = ii;
			sign_pv		= 0;
			datai.resize((((nstop-1)*_fsh10+_flen-_fsh10)-((nstart-1)*_fsh10+1))+1);
			datai = dfdata.segment((nstart-1)*_fsh10,(((nstop-1)*_fsh10+_flen-_fsh10)-((nstart-1)*_fsh10+1))+1); 
			
			for(int jj=nstart;jj <= nstop-1;jj++)
			{
				for(int h=1; h<=_flen; h++)
				{
					e[jj-1] = e[jj-1]+datai[(jj-nstart)*_fsh10+h-1]*datai[(jj-nstart)*_fsh10+h-1];
				}
				if(e[jj-1] <= ENERGYFLOOR)
					e[jj-1] = ENERGYFLOOR;
			}
			e[nstop - 1] = e[nstop-2];
			e_temp.resize(nstop-nstart+1); eY.resize(nstop-nstart+1);
			e_temp = e.segment(nstart-1,nstop-nstart+1);
			sort_up_np(e_temp,eY);
			emin = eY[int((nstop-nstart+ 1)* 0.1)- 1];

			for(int jj=nstart + 1;jj <= nstop; jj++)
			{
				postsnr[jj - 1] = log10(e[jj-1])-log10(emin);
				if(postsnr[jj-1] < 0)
				   postsnr[jj - 1] = 0;
				D[jj-1] = sqrt(abs(e[jj-1]-e[jj-2])*postsnr[jj-1]);
			}
			D[nstart-1]	  = D[nstart];
			int Dexp_num  = Dexpl+(nstop-nstart+1)+Dexpr;
			Dexp.resize(Dexp_num);
			Dexp << D[nstart-1]*VectorXd::Ones(Dexpl),D.segment(nstart-1,nstop-nstart+1),D[nstop-1]*VectorXd::Ones(Dexpr);
			
			for(int jj=0; jj <=(nstop-nstart);jj++)
			{
				Dexp_temp = Dexp.segment(jj,Dexpl+Dexpr);
				Dsmth[nstart+jj-1] = Dexp_temp.sum();
			}

			double Dsmth_thres = (Dsmth.segment(nstart-1,nstop-nstart+1).array()*
			pv01.segment(nstart-1,nstop-nstart+1).array()).sum()/(pv01.segment(nstart-1,nstop-nstart+1).sum());
			for(int jj = nstart;jj <= nstop; jj++)
			{
				if(Dsmth[jj - 1] > Dsmth_thres * vadThres * snr_adj)
					snre_vad[jj-1] = 1;
			}
		}
	}
	    pv_vad      = snre_vad;   
	int	nexpl		= 33;
	int	nexpr		= 47;
	int	sign_vad	= 0;
	for(int ii=1; ii<= _nfr10; ii++)
	{
		if((snre_vad[ii - 1] == 1) && (sign_vad == 0))
		{
			nstart		= ii;
			sign_vad	= 1;
		}
		else if((snre_vad[ii-1] == 0 || ii == _nfr10) && (sign_vad == 1))
		{
			nstop = ii - 1;
			if(ii == _nfr10)  nstop = ii;
			sign_vad = 0;
			int j = nstart;
			for(j = nstart; j <= nstop; j++)
			{
				if(pv01[j - 1] == 1)
					break;
			}
			//pv_vad.segment(nstart-1, max(j-nexpl-1,1)-nstart+1).fill(0);
			for(int nn = nstart; nn <= max(j-nexpl-1,1); nn++)
				pv_vad[nn - 1] = 0;
			j = 0;			
			for(j=0;j <= (nstop-nstart);j++)
			{
				if(pv01[nstop-j-1] == 1)
					break;
			}
			for(int nn = nstop-j+1+ nexpr; nn <= nstop; nn++ )
				pv_vad[nn-1] = 0;
			j = 0;
		}
	}	
	/* the second round */
	nexpl		= 5;
	nexpr		= 12;
	sign_vad	= 0;
	for(int ii=1; ii <= _nfr10; ii++)
	{
		if((snre_vad[ii - 1] == 1) && (sign_vad == 0))
		{
			nstart		= ii;
			sign_vad	= 1;
		}
		else if((snre_vad[ii - 1] == 0 || ii == _nfr10) && (sign_vad == 1))
		{
			nstop = ii - 1;
			if ( ii == _nfr10 ) nstop = ii;
			sign_vad	= 0;

			if( pv01.segment(nstart-1,nstop-nstart+1).sum() >4)
			{
				int j = nstart;
				for (j = nstart; j <= nstop; j++)
				{
					if(pv01[j-1] == 1)  break;
				}

				for(int nn = max(j-nexpl,1); nn <= j-1; nn++ )
					pv_vad[nn - 1] = 1;
				j = 0;
				for(j = 0;j <= nstop - nstart;j++ )
				{
					if(pv01[nstop - j - 1] == 1)
						break;
				}
				for (int nn = nstop - j + 1; nn <= min( nstop - j + nexpr, _nfr10 ); nn++)
					pv_vad[nn - 1] = 1;
			}			
			double esegment = e.segment(nstart-1,nstop-nstart+1).sum()/(nstop-nstart+1);
			if(esegment < 0.001)
			{
				pv_vad.segment(nstart-1,nstop-nstart+1).fill(0);
			}
			if( pv01.segment(nstart-1,nstop-nstart+1).sum() <= 2)
			{
				pv_vad.segment(nstart-1,nstop-nstart+1).fill(0);
			}
		}
	}
	sign_vad    = 0;
	double esum = 0;
	for(int ii = 1;ii <= _nfr10; ii++)
	{
		if(pv_vad[ii - 1] == 1 && sign_vad == 0)
		{
			nstart		= ii;
			sign_vad	= 1;
		}
		else if((pv_vad[ii - 1] == 0 || ii == _nfr10) && sign_vad == 1)
		{
			nstop = ii - 1;
			if(ii == _nfr10) nstop = ii;
			sign_vad	= 0;
			esum = esum + e.segment(nstart-1,nstop-nstart+1).sum();
		}
	}

	double eave = esum/pv_vad.sum();	
	sign_vad = 0;
	for(int ii = 1;ii <= _nfr10; ii++)
	{
		if(pv_vad[ii - 1] == 1 && sign_vad == 0)
		{
			nstart		= ii;
			sign_vad	= 1;
		}
		else if((pv_vad[ii - 1] == 0 || ii == _nfr10) && (sign_vad == 1))
		{
			nstop = ii - 1;
			if(ii == _nfr10)  nstop = ii;
			sign_vad	= 0;
			if(e.segment(nstart-1,nstop-nstart+1).sum()/(nstop-nstart+1) < eave*0.05)
			{
				pv_vad.segment(nstart-1,nstop-nstart+1).fill(0);
			}
		}
	}
	int ii = 0 ;
}
void VAD::vad_snre_round1(VectorXd& dfdata,VectorXi& snre_vad,VectorXd& e)
{
	int		 Dexpl    = 10;
	int		 Dexpr    = 10;
	double   vadThres = 0.25;
	e        = VectorXd::Zero(_nfr10);
	for(int i=1;i<=_nfr10;i++)
	{
		for(int j=1;j<=_flen;j++)
		   e(i-1) = e(i-1)+dfdata((i-1)*_fsh10+j-1)*dfdata((i-1)*_fsh10+j-1);
		if(e(i-1)<ENERGYFLOOR)
		   e(i-1) = ENERGYFLOOR;
	}
	VectorXd emin     = VectorXd::Ones(_nfr10);
	int NESEG = 200;
	if(_nfr10 < NESEG) NESEG = _nfr10;
	VectorXd eY;
	int i;
	for(i=1; i<=floor(double(_nfr10)/NESEG);i++)
	{
		sort_up_np(e.segment((i-1)*NESEG,NESEG),eY);
		emin.segment((i-1)*NESEG,NESEG).fill(eY(floor(NESEG*0.1)-1)); 
		if( i != 1)
		   emin.segment((i-1)*NESEG,NESEG).fill(0.9*emin((i-1)*NESEG-1)+0.1*emin((i-1)*NESEG));    
	}
	eY.resize(0);
    i = floor(double(_nfr10)/NESEG);
	if (i*NESEG != _nfr10)
	{
		sort_up_np(e.segment((i-1)*NESEG,_nfr10-(i-1)*NESEG),eY);
	    emin.segment(i*NESEG, _nfr10-i*NESEG).fill(eY(floor((_nfr10-(i-1)*NESEG)*0.1)-1));  
	    emin.segment(i*NESEG,_nfr10-i*NESEG).fill(0.9*emin(i*NESEG-1)+0.1*emin(i*NESEG)); 
	}
	VectorXd D       =  VectorXd::Zero(_nfr10);
	VectorXd postsnr =  VectorXd::Zero(_nfr10);
	for(i = 2;i<= _nfr10;i++)
	{
		postsnr(i-1) =log10(e(i-1))-log10(emin(i-1));
		if(postsnr(i-1)<0)   postsnr(i-1) = 0;
	    D(i-1)=sqrt(abs(e(i-1)-e(i-2))*postsnr(i-1)); 
	}
	D(0) = D(1);
	VectorXd  Dexp(Dexpl+D.size()+Dexpr);
	Dexp << D(0)*VectorXd::Ones(Dexpl).array(),D,D(_nfr10-1)*VectorXd::Ones(Dexpr).array();
	VectorXd  Dsmth      = VectorXd::Zero(_nfr10);
	for( i = 1;i<=_nfr10;i++)
		Dsmth(i-1) = Dexp.segment(i-1,Dexpl+Dexpr+1).sum();
	VectorXd  Dsmth_max	 = VectorXd::Zero(_nfr10);
	
	for( i = 1;i<=floor(double(_nfr10)/NESEG);i++)
		Dsmth_max.segment((i-1)*NESEG,NESEG).fill(Dsmth.segment((i-1)*NESEG,NESEG).maxCoeff());
	i = floor(double(_nfr10)/NESEG);
	if(i*NESEG != _nfr10)
		Dsmth_max.segment(i*NESEG,_nfr10-i*NESEG).fill(Dsmth.segment((i-1)*NESEG,_nfr10-(i-1)*NESEG).maxCoeff());
	//VectorXi 
	snre_vad = VectorXi::Zero(_nfr10);
	for( i = 1;i<=_nfr10;i++)
	{
		if( Dsmth(i-1) > Dsmth_max(i-1) * vadThres)
			snre_vad(i-1) = 1;
	}

}
void VAD::pitchblockdetect(VectorXd& pv01,VectorXd& pitch,VectorXd& pvblk)
{
	if(_nfr10 == pv01.size()+1)
		pv01(_nfr10-1) = pv01(_nfr10-2);
	int sign_pv = 0;
	int nstart = 0,nstop = 0;
	VectorXd pitchseg;

	for(int i=1;i<=_nfr10;i++)
	{
		if((pv01(i-1)==1) && (sign_pv ==0))
		{
			nstart  = i;
			sign_pv = 1; 
		}
		else if((pv01(i-1)==0 || i ==_nfr10)&&(sign_pv == 1))
		{
			nstop = i;
		    if( i==_nfr10)
				nstop = i+1;
			sign_pv = 0;
			pitchseg.resize(nstop-nstart);
			pitchseg = VectorXd::Zero(nstop-nstart);
		    for(int j=nstart;j<=nstop-1;j++)
				pitchseg(j-nstart) = pitch(j-1);
			pitchseg = pitchseg - pitchseg.mean()*VectorXd::Ones(pitchseg.size());
			for(int jj=0;jj<pitchseg.size();jj++)
				pitchseg(jj) = abs(int(pitchseg(jj)+0.5));
			if((pitchseg.sum()==0) && (nstop-nstart+1>=10))
				pv01.segment(nstart-1,nstop-nstart).fill(0);
		}
	}
	sign_pv  = 0;
	pvblk = pv01;
	for(int i=1;i<= _nfr10;i++)
	{
		if((pv01(i-1)==1)&&(sign_pv == 0))
		{
			nstart  = i;
		    sign_pv = 1;
			pvblk.segment(max(nstart-60,1)-1, nstart-max(nstart-60,1)+1).fill(1);
		}
		else if((pv01(i-1)==0 || i ==_nfr10)&&(sign_pv==1))
		{
			nstop   = i;
			sign_pv = 0;
		    pvblk.segment(nstop-1,min(nstop+60,_nfr10)-nstop+1) .fill(1);
		}
	}
}
/* VAD smooth */
void VAD::merge(double* &vseg, int vseg_num, double* &vseg_mer, int &vseg_mer_num)   /* vseg(vseg_num,2) */
{
	double* vseg_merge = new double[vseg_num * 2]();
	for(int ii = 0;ii<vseg_num;ii++)
	{
		vseg_merge[2*ii+0] = vseg[2*ii+ 0];
		vseg_merge[2*ii+1] = vseg[2*ii+ 1];
	}
	for ( int ii = 2; ii <= vseg_num; ii++ )
	{
		if ( (vseg[2 * (ii - 1) + 0] - vseg[2 * (ii - 1 - 1) + 1]) <= 0.5 )
		{
			vseg_merge[2 * (ii - 1) + 0]	= 0;
			vseg_merge[2 * (ii - 2) + 1]	= 0;
		}
	}
	int t1 = 1, t2 = 1;
	for ( int ii = 1; ii <= vseg_num; ii++ )
	{
		if ( vseg_merge[2 * (ii - 1) + 0] != 0 )
			t1++;
		if ( vseg_merge[2 * (ii - 1) + 1] != 0 )
			t2++;
	}
	vseg_mer_num	= max( t1 - 1, t2 - 1 );
	vseg_mer	= new double[vseg_mer_num * 2]();
	t1		= 1; t2 = 1;
	for ( int ii = 1; ii <= vseg_num; ii++ )
	{
		if ( vseg_merge[2 * (ii - 1) + 0] != 0 )
		{
			vseg_mer[2 * (t1 - 1) + 0] = vseg_merge[2 * (ii - 1) + 0];
			t1++;
		}
		if ( vseg_merge[2 * (ii - 1) + 1] != 0 )
		{
			vseg_mer[2 * (t2 - 1) + 1] = vseg_merge[2 * (ii - 1) + 1];
			t2++;
		}
	}
	delete[]vseg_merge; vseg_merge = NULL;
}
/* VAD split */
void VAD::split( double* &vseg_mer, int vseg_mer_num, double* &vseg_sp)   /* vseg_mer(vseg_mer_num,2) */
{
	vseg_sp = new double[vseg_mer_num * 2]();
	for ( int ii = 0; ii < vseg_mer_num; ii++ )
	{
		vseg_sp[2 * ii + 0]	= vseg_mer[2 * ii + 0];
		vseg_sp[2 * ii + 1]	= vseg_mer[2 * ii + 1];
	}
	for ( int ii = 2; ii <= vseg_mer_num; ii++ )
	{
		if ((vseg_sp[2 * (ii - 1) + 0] - vseg_sp[2 * (ii - 1 - 1) + 1]) <= 0.5 )
			cout << " wrong in merge" << endl;
	}
	double gap = 0;
	for ( int ii = 2; ii <= vseg_mer_num; ii++ )
	{
		if ( (vseg_sp[2 * (ii - 1) + 0] - vseg_sp[2 * (ii - 1 - 1) + 1]) <= 1 )
		{
			gap				= vseg_sp[2 * (ii - 1) + 0] - vseg_sp[2 * (ii - 1 - 1) + 1];
			vseg_sp[2 * (ii - 1) + 0]	= vseg_sp[2 * (ii - 1) + 0] - gap / 2.0;
			vseg_sp[2 * (ii - 1 - 1) + 1]	= vseg_sp[2 * (ii - 1 - 1) + 1] + gap / 2.0;
		}
		else
		{
			vseg_sp[2 * (ii - 1) + 0]	= vseg_sp[2 * (ii - 1) + 0] - 0.2;
			vseg_sp[2 * (ii - 1 - 1) + 1]	= vseg_sp[2 * (ii - 1 - 1) + 1] + 0.2;
		}
	}
	vseg_sp[0] = vseg_sp[0] - 0.2;
	if ( vseg_sp[0] < 0 )
		vseg_sp[0] = 0;
	vseg_sp[2 * (vseg_mer_num - 1) + 1] = vseg_sp[2 * (vseg_mer_num - 1) + 1] + 0.2;
}


