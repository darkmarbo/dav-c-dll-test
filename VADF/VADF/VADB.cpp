#include "VAD.h"
#include "VADB.h"

double* VADB(short* buf,double SampleRate,int SamplePoints,int& seg_num)
{
	//double *vad_seg = NULL;                                    //��ʱ���������� 
	//VAD EPitch(buf, SamplePoints, SampleRate, vad_seg, seg_num);   //����vad_seg������
    //return vad_seg;


	////////////////////////////////////////// �޸ĺ� 
	int time_small = 60; // ÿ�δ���Ĳ������� ����Ӧ2���ӵ����ݣ�

	int smp_rate_int = int(SampleRate);  // ������
	int smp_count = SamplePoints;  // ��ȡ�����ļ��ܳ��� ������֡��
	int num = 0;
	double *seg_res = new double[2000];
	double *seg_res_now = seg_res;
	short *buff_tmp = buf;  // ָ������buff���ݵ�ָ��
	seg_num = 0;

	int SMP_CNT_SMALL = smp_rate_int * time_small;  // ÿ�δ���Ĳ������� 2���� 

	int smp_cnt_small = 0;  // ʵ���ϴ���Ĳ������� 
	double *seg_small = new double[1000]; // �洢ÿ��С���зֶ����� ���ص��зֽ��  
	
	int smp_count_remain = smp_count;  // ʣ�� δ�зֵ����ݶΣ������㣩
	int flag_over = 0;  // �Ƿ����������������� 
	double time_frn = 0;  // ��ǰ��ʱ���Ƕ��� 
	int smp_cnt_frn = 0;
	int num_frn = 0;


	// ѭ������ÿ��С�� 
	while (smp_count_remain > 0)
	{
		short *buff_small = NULL;   // �洢ÿ��Ҫ�����buff

		if (smp_count_remain > SMP_CNT_SMALL)
		{
			smp_cnt_small = SMP_CNT_SMALL;
		}
		else
		{
			flag_over = 1;
			smp_cnt_small = smp_count_remain;
		}

		buff_small = new short[smp_cnt_small];
		memcpy(buff_small, buff_tmp, sizeof(short)*smp_cnt_small);

		VAD EPitch(buff_small,  smp_cnt_small, SampleRate, seg_small, num);


		// �Ȱ�ǰ num-1�������� 
		for (int ii = 0; ii < num - 1; ii++)
		{
			seg_small[2 * ii] += time_frn;
			seg_small[2 * ii + 1] += time_frn;
			seg_num++;

		}		
		memcpy(seg_res_now, seg_small, sizeof(double)*(2 * (num - 1)));
		seg_res_now += (2 * (num - 1));  // �洢��ǰ��num - 1�����׵�

		// �ж��Ƿ������һ��  
		//printf("seg_num %d\n", num);
		if (flag_over == 1)
		{
			seg_num ++;
			seg_small[2 * (num - 1)] += time_frn;
			seg_small[2 * (num - 1) + 1] += time_frn;

			memcpy(seg_res_now, seg_small+2*(num-1), sizeof(double)* 2);
			seg_res_now += 1;
			break;
		}
 

		// ����ʣ��λ��ж೤ �����￪ʼ
		
		time_frn = seg_small[2 * (num - 1) - 1]; // ǰ�洦����� ���е�ʱ�䳤�� 
		smp_cnt_frn = smp_rate_int*time_frn; // ǰ�洦����� ���еĲ�����

		smp_count_remain = smp_count - smp_cnt_frn;
		buff_tmp = buf + smp_cnt_frn;

	}

	// �ͷ� �� �ر� 
	return seg_res;

}


