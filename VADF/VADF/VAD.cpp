#include "VAD.h"

/* Class Construct Function */
VAD::VAD(short* buf,int length,double sampleRate,double* &vad_seg_dy,int &vad_num_dy)
{
 	/* data transform */
	_fs         = sampleRate;					//采样率
	_nSamples   = length;						//采样长度
	_data.resize(_nSamples);
	for(int i = 0; i < _nSamples; i++ )
		_data[i] = (double) (buf[i])/32768;     //数据转换成功
	_flen	= (int) _fs * 0.025;
	_fsh10	= (int) _fs * 0.01;
	_nfr10	= (int) floor( double( (_nSamples - (_flen - _fsh10) ) / _fsh10));

	VectorXd  pv01,pitch;
	pitchestm(_data,pv01,pitch);	
	double    b[2] = {0.9770,-0.9770};
	double    a[2] = {1.000,-0.9540};
	VectorXd  fdata = VectorXd::Zero(_data.size());
	fdata(0)  = b[0] * _data(0)/a[0];
	for ( int ii = 1; ii < _nSamples; ii++ )
		fdata(ii) = (b[0]*_data(ii)+b[1]*_data(ii-1)-a[1]*fdata(ii-1))/a[0];
	_data.resize(0);
	VectorXd pvblk;
	pitchblockdetect(pv01,pitch,pvblk);
	VectorXi snre_vad1;  VectorXd e;
	vad_snre_round1(fdata,snre_vad1,e);
	
	int      sign_vad      = 0; 
	int      nstart        = 0;
	int      nstop         = 0;
	int      n_noise_samp  = 0;
	VectorXd noise_seg       =  VectorXd::Zero(floor(double(_nfr10)/1.6));
	MatrixXd noise_samp_temp =  MatrixXd::Zero(_nfr10,2);
	for(int i=1; i<=_nfr10;i++)
	{
		if(snre_vad1(i-1)==1  && sign_vad==0)
		{
			sign_vad = 1;
			nstart   = i;
		}
		else if((snre_vad1(i-1)==0 || i==_nfr10) && (sign_vad == 1))
		{
			sign_vad = 0;
			nstop    = i-1;
			if(pv01.segment(nstart-1,nstop-nstart+1).sum() == 0)
			{
				noise_seg.segment(int(nstart/1.6+0.5)-1, floor(nstop/1.6)-int(nstart/1.6+0.5)+1).fill(1);
				++n_noise_samp;
				noise_samp_temp(n_noise_samp-1,0) = (nstart-1)*_fsh10+1;
				noise_samp_temp(n_noise_samp-1,1) =  nstop * _fsh10;
			}
		}
	}
	MatrixXd noise_samp = noise_samp_temp.block(0,0,n_noise_samp,2);
	noise_samp_temp.resize(0,0);
	VectorXd dfdatarm;
    specsub_rats_noiseseg_1fn(fdata,noise_seg,pv01,dfdatarm);
	for(int ii=1;ii<= n_noise_samp;ii++)
	{
	   dfdatarm.segment(noise_samp(ii-1,0)-1,noise_samp(ii-1,1)- noise_samp(ii-1,0)+1).fill(0);
	}
	double  vadThres = 0.1;
	VectorXd pv_vad2;
	vad_snre_pv(dfdatarm,pv01,pvblk,vadThres,pv_vad2);
	/* -------------------------------------------------------- */
	/* 计算最终结果 vad_seg的大小*/
	sign_vad = 0;
	int n_vad_seg = 0;
	for(int ii=1;ii<= _nfr10;ii++)
	{
		if((pv_vad2[ii - 1] == 1) && (sign_vad == 0))
		{
			nstart = ii; sign_vad = 1;
		}
		else if((pv_vad2[ii - 1] == 0 || ii == _nfr10) && sign_vad == 1)
		{
			nstop = ii - 1; sign_vad = 0; n_vad_seg = n_vad_seg + 1;
		}
	}
	MatrixXd vad_seg(n_vad_seg,2);     //the end data save
    n_vad_seg = 0;
	sign_vad  = 0;
	nstart	  = nstop = 0;
	for( int ii = 1; ii <= _nfr10; ii++)
	{
		if((pv_vad2[ii-1] == 1) && (sign_vad == 0))
		{
			nstart		= ii;
			sign_vad	= 1;
		}
		else if((pv_vad2[ii-1] == 0 || ii == _nfr10) && sign_vad == 1)
		{
			nstop					  = ii - 1;
			sign_vad				  = 0;
			n_vad_seg				  = n_vad_seg + 1;
			vad_seg(n_vad_seg - 1,0)  = nstart * 0.01;
			vad_seg(n_vad_seg - 1,1)  = nstop  * 0.01;
		}
	}
	/* convert and save the vad outcome */
	vad_num_dy = vad_seg.rows();
    vad_seg_dy = new double[vad_seg.rows()*vad_seg.cols()]();
	for(int ii=0;ii<vad_seg.rows();ii++)
	{
		vad_seg_dy[2 * ii + 0]	= vad_seg(ii,0);
		vad_seg_dy[2 * ii + 1]	= vad_seg(ii,1);
	}
}
