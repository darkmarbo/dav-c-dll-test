#include "VAD.h"
#include "VADB.h"

double* VADB(short* buf,double SampleRate,int SamplePoints,int& seg_num)
{
	//double *vad_seg = NULL;                                    //暂时保存在这里 
	//VAD EPitch(buf, SamplePoints, SampleRate, vad_seg, seg_num);   //返回vad_seg储存结果
    //return vad_seg;


	////////////////////////////////////////// 修改后 
	int time_small = 60; // 每次处理的采样点数 （对应2分钟的数据）

	int smp_rate_int = int(SampleRate);  // 采样率
	int smp_count = SamplePoints;  // 读取语音文件总长度 （多少帧）
	int num = 0;
	double *seg_res = new double[2000];
	double *seg_res_now = seg_res;
	short *buff_tmp = buf;  // 指向整个buff数据的指针
	seg_num = 0;

	int SMP_CNT_SMALL = smp_rate_int * time_small;  // 每次处理的采样点数 2分钟 

	int smp_cnt_small = 0;  // 实际上处理的采样点数 
	double *seg_small = new double[1000]; // 存储每个小的切分段数据 返回的切分结果  
	
	int smp_count_remain = smp_count;  // 剩余 未切分的数据段（采样点）
	int flag_over = 0;  // 是否处理完了整个大语音 
	double time_frn = 0;  // 当前的时间是多少 
	int smp_cnt_frn = 0;
	int num_frn = 0;


	// 循环处理每个小段 
	while (smp_count_remain > 0)
	{
		short *buff_small = NULL;   // 存储每次要处理的buff

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


		// 先把前 num-1个拷贝了 
		for (int ii = 0; ii < num - 1; ii++)
		{
			seg_small[2 * ii] += time_frn;
			seg_small[2 * ii + 1] += time_frn;
			seg_num++;

		}		
		memcpy(seg_res_now, seg_small, sizeof(double)*(2 * (num - 1)));
		seg_res_now += (2 * (num - 1));  // 存储了前面num - 1个靠谱的

		// 判断是否是最后一个  
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
 

		// 计算剩余段还有多长 从哪里开始
		
		time_frn = seg_small[2 * (num - 1) - 1]; // 前面处理过的 所有的时间长度 
		smp_cnt_frn = smp_rate_int*time_frn; // 前面处理过的 所有的采样点

		smp_count_remain = smp_count - smp_cnt_frn;
		buff_tmp = buf + smp_cnt_frn;

	}

	// 释放 与 关闭 
	return seg_res;

}


