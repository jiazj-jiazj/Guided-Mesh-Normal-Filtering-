
#include "ExKernelT.h"
#include <iostream>
#include <iterator>
#include <random>

namespace MeshN { 

	///////////////////////////////////////////////////////////////////////////////
	// Implementation of member functions of ExKernelT
	/////////////////////////////////////////////////////////////////////////////// 
	template<class ExItems>
	ExKernelT<ExItems>::ExKernelT(): KernelT<ExItems>() {

		//kdtree_ = NULL;
		//ps_     = NULL;
		isNormal_ = false;
		isArea_ = false;

	}
	//////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	ExKernelT<ExItems>::~ExKernelT(){

		//if(kdtree_ != NULL) delete kdtree_;

		//if(ps_ != NULL)     delete ps_;

	}
	//////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	void ExKernelT<ExItems>::meshInit(){

		update_facet_normals();//��������������Ƭ�ķ���
		update_area();
		update_edge_length();//������������ı߳�


	}

	///////////////////////////////////////////////////////////////////////////////
	template <class ExItems>
	double ExKernelT<ExItems>::count_K(Normal normal_1, Normal normal_2, double egama) {    //����˺���ֵ����������������һ���˲���

		Normal normal_3 = normal_1 - normal_2;
		double length = sqrt(normal_3.data_[0] * normal_3.data_[0] + normal_3.data_[1] * normal_3.data_[1] + normal_3.data_[2] * normal_3.data_[2]);
		double k = exp(-(length*length) / (2 * (egama*egama)));
		return k;
	};
	template<class ExItems>
	void
		ExKernelT<ExItems>::adjustFaceNormal_2015() {    //ʵ��2015������  Guided Mesh Normal Filtering

		double sigama_r;  //�˺�������sigama_r
		std::cout << "������sigama_r [0.2-0.6] ";
		std::cin >> sigama_r;
		double sigama_s = 0;
		int number = 0;
		FacetIterator fit = facet_begin();
		Coord sum_centroid(0, 0, 0);
		Normal sum_normal(0, 0, 0);
		for (; fit != facet_end(); fit++) {
			number++;
			std::vector<FacetHandle> fhs;
			HalfedgeHandle hh = fit->halfedge_handle_;
			FacetHandle    fh = facet_handle(hh);

			Coord sum_centroid(0, 0, 0);
			Normal sum_normal(0, 0, 0);

			Coord coord = calc_centroid(fh);
			Normal normal = calc_normal(fh);

			sum_centroid += coord;
			sum_normal += normal;
		}  //�����׼��
		Coord aver_centroid = sum_centroid / number;   //��Ƭ��ƽ������
		Normal aver_normal = sum_normal / number;     //��Ƭ��ƽ������

		fit = facet_begin();
		for (; fit != facet_end(); fit++) {    //���sigama_s
			number++;
			std::vector<FacetHandle> fhs;
			HalfedgeHandle hh = fit->halfedge_handle_;
			FacetHandle    fh = facet_handle(hh);

			Coord coord = calc_centroid(fh);
			Normal normal = calc_normal(fh);

			Coord normal_3 = coord - aver_centroid;
			double length_s = (normal_3.data_[0] * normal_3.data_[0] + normal_3.data_[1] * normal_3.data_[1] + normal_3.data_[2] * normal_3.data_[2]);
			sigama_s += length_s;
		}
		sigama_s = sqrt(sigama_s / (number - 1));

		int facet_num = facet_size();

		int iteration;
		std::cout << "Input facenormal filtering numbers (10-75��): ";
		std::cin >> iteration;
		std::cout << "���ڵ������ȽϺ�ʱ��" << std::endl;
		std::cout << "Please wait.....  " << std::endl;
		int i = 0, j = 0;
		clock_t t1 = clock();
		do {                        //��ʼ����

			std::vector<Normal> updateFacetNormal;   //������Ƭ����ķ�����
			updateFacetNormal.resize(facet_num);

			std::vector<Normal> updateFacet_n_;   //����patch�鷨��ı�ѡ������
			updateFacet_n_.resize(facet_num);

			std::vector<Normal> updateFacet_g_;   //������Ƭ��������������
			updateFacet_g_.resize(facet_num);

			double * h_min = (double*)malloc(sizeof(double) * 200000);
			memset(h_min, 127, sizeof(h_min));

			double *sita = (double*)malloc(sizeof(double) * 200000); //���ÿ��patch�Ħ���P��
			double * h = (double*)malloc(sizeof(double) * 200000);   //���ÿ��patch��H��P��=����P��*R��P��
			double * r = (double*)malloc(sizeof(double) * 200000);    //���ÿ��patch��R��P��

			memset(sita, 0, sizeof(sita));
			memset(h, 0, sizeof(h));
			memset(r, 0, sizeof(r));

			FacetIterator fit = facet_begin();
			for (; fit != facet_end(); fit++) {          //���� ���ÿ��patch�����ֵ
				std::vector<FacetHandle> fhs;
				HalfedgeHandle hh = fit->halfedge_handle_;
				FacetHandle    fh = facet_handle(hh);
				fhs.push_back(fh);

				std::vector<VertexHandle> varry;
				varry.resize(3);
				varry[0] = vertex_handle(hh);   //�õ�һ�������������
				varry[1] = vertex_handle(prev_halfedge_handle(hh));
				varry[2] = vertex_handle(next_halfedge_handle(hh));


				for (i = 0; i < 3; i++) {        //����Ƭ��Χ�Ĺ�������Ƭ����patch����
					HalfedgeHandle vahh = halfedge_handle(varry[i]);
					HalfedgeHandle vacss(vahh);
					do {
						FacetHandle vafh = facet_handle(vacss);
						if (vafh.is_valid()) {

							for (j = 0; j < fhs.size(); j++)
								if (vafh == fhs[j]) break;
							if (j >= fhs.size()) fhs.push_back(vafh);
						}

						vacss = cw_rotated(vacss);

					} while (vacss != vahh);
				}

				double sum = 0;

				Normal normal_t(0, 0, 0); //hatch��ѡ��������

				double *together = (double*)malloc(sizeof(double) * 400000);

				double max_saliency = 0;
				double sum_saliency = 0;

				for (int i = 0; i < fhs.size(); i++)   //��������ĳһ��Ƭ��������Ƭ
				{
					int idx = fhs[i].idx();
					float area = calc_facet_area(fhs[i]);
					Normal norr = calc_normal(fhs[i]);
					normal_t += (area*norr);                   //��ѡ�����������������Ȩ��

					std::vector<FacetHandle> Face_common_edge;  //����ÿ����Ƭ���ߵ���Ƭ

					getNeighborFaceN1(fhs[i], Face_common_edge);  //Face_common_edge ��ȡ��Ƭ���ߵ�������Ƭ

					for (int k = 0; k < Face_common_edge.size(); k++)  //����ÿ����Ƭ�������Ƭ
					{

						if (!together[(fhs[i].idx() + Face_common_edge[k].idx()) ])    //����ù����û�б����ʹ�
						{
							Normal normal_together_1 = calc_normal(fhs[i]);
							Normal normal_together_2 = calc_normal(Face_common_edge[k]);

							Normal normal_together_3 = normal_together_1 - normal_together_2;
							double length_s = sqrt(normal_together_3.data_[0] * normal_together_3.data_[0] + normal_together_3.data_[1] * normal_together_3.data_[1] + normal_together_3.data_[2] * normal_together_3.data_[2]);

							sum_saliency += length_s;

							if (max_saliency < length_s)
							{
								max_saliency = length_s;   //ÿ�θ���max_saliency 
							}

							together[fhs[i].idx() + Face_common_edge[k].idx()] = 1;   //��¼�ù�����Ѿ������ʹ����´β��ٷ���
						}
					}

					for (int j = 0; j < fhs.size(); j++)
					{
						Normal normal_1 = calc_normal(fhs[i]);
						Normal normal_2 = calc_normal(fhs[j]);
						Normal normal_3 = normal_1 - normal_2;
						double length = sqrt(normal_3.data_[0] * normal_3.data_[0] + normal_3.data_[1] * normal_3.data_[1] + normal_3.data_[2] * normal_3.data_[2]);
						sum += length;
						if (length > sita[fh.idx()])
							sita[fh.idx()] = length;   //ÿ�θ��¸���Ƭ�Ħ�(P)
					}

				}
				free(together);
				if (sum_saliency)
					r[fh.idx()] = max_saliency / sum_saliency;  //ÿ�θ��¸���Ƭ��R(P)

				updateFacet_n_[fh.idx()] = normal_t.normalize();     //ÿ�θ��¸���Ƭ��ƽ������

				h[fh.idx()] = sita[fh.idx()] * r[fh.idx()];   //���������Ƭ��H(P)ֵ
				for (int i = 0; i < fhs.size(); i++)
				{
					FacetHandle facehandle = fhs[i];
					if (h[fh.idx()] < h_min[facehandle.idx()])     //����h(p)����ø���Ƭ����
					{
						h_min[facehandle.idx()] = h[fh.idx()];
						updateFacet_g_[facehandle.idx()] = updateFacet_n_[fh.idx()];//ÿ�θ���ÿ����Ƭ��������������
					}
				}
			}

			fit = facet_begin();
			for (; fit != facet_end(); fit++) {       //����ÿ��patch�飬����ÿ����Ƭ�����������������Ƭ�ĸ��·�����
				std::vector<FacetHandle> fhs;
				HalfedgeHandle hh = fit->halfedge_handle_;
				FacetHandle    fh = facet_handle(hh);

				Normal centroid_i = calc_centroid(fh);
				Normal normal_g_i = updateFacet_g_[fh.idx()];

				fhs.push_back(fh);
				//VertexHandle varry[3];
				std::vector<VertexHandle> varry;
				varry.resize(3);
				varry[0] = vertex_handle(hh);
				varry[1] = vertex_handle(prev_halfedge_handle(hh));
				varry[2] = vertex_handle(next_halfedge_handle(hh));
				for (i = 0; i < 3; i++) {                        
					HalfedgeHandle vahh = halfedge_handle(varry[i]);
					HalfedgeHandle vacss(vahh);
					do {
						FacetHandle vafh = facet_handle(vacss);
						if (vafh.is_valid()) {

							for (j = 0; j < fhs.size(); j++)
								if (vafh == fhs[j]) break;
							if (j >= fhs.size()) fhs.push_back(vafh);
						}

						vacss = cw_rotated(vacss);

					} while (vacss != vahh);
				}

				Normal normal_face(0, 0, 0);
				for (int j = 0; j < fhs.size(); j++)
				{
					Normal normal_j = calc_normal(fhs[j]);
					Normal centroid_j = calc_centroid(fhs[j]);
					Normal normal_g_j = updateFacet_g_[fhs[j].idx()];
					float s = calc_facet_area(fhs[j]);                     //������Ƭ�����
					float ks = count_K(centroid_i, centroid_j, sigama_s);  //����˺�����ֵ
					float kr = count_K(normal_g_i, normal_g_j, sigama_r);   //����˺�����ֵ
					normal_face += s * ks*kr*normal_j;
				}

				updateFacetNormal[fh.idx()] = normal_face.normalize();  //������׼��
			}

			for (fit = facet_begin(); fit != facet_end(); fit++) {                     //������Ƭ������
				Normal normal = updateFacetNormal[facet_handle(fit->halfedge_handle_).idx()];

				facet_ref(facet_handle(fit->halfedge_handle_)).normal_ = updateFacetNormal[facet_handle(fit->halfedge_handle_).idx()].normalize();

				Normal normall = facet_ref(facet_handle(fit->halfedge_handle_)).normal_;
			}

			free(h_min);   //�ͷŶ��Ͽ��Ŀռ�
			free(sita);
			free(h);
			free(r);
		} while (--iteration);

		clock_t t2 = clock();
		std::cout << "The time of Sun's normal filter: " << (t2 - t1)*1.0 / CLOCKS_PER_SEC << "s" << std::endl;
	}


	template<class ExItems>
	void ExKernelT<ExItems>::Guided_Mesh_Normal_Filtering_2015(){

		adjustFaceNormal_2015();    //������Ƭ����
		//adjustFaceNormal_Based_on_FE();
		int vertex_num = vertex_size();
		int iterations;
		std::cout << "Input vertex update iterations (10-30��): ";
		std::cin >> iterations;
		std::cout << "���ڵ������ȽϺ�ʱ��" << std::endl;
		std::cout << "Please wait.....  " << std::endl;
		int i = 0;

		clock_t t1 = clock();
		do {
			std::vector<Coord> updateVertexPosition;
			updateVertexPosition.resize(vertex_num);
			for (i = 0; i < vertex_num; i++) {
				VertexHandle vh(i);
				Coord&       vc = coord(vh);
				HalfedgeHandle& hh = halfedge_handle(vh);
				HalfedgeHandle  css(hh);
				do {
					HalfedgeHandle opp_hh = opposite_halfedge_handle(css);
					Coord&         opp_vc = coord(vertex_handle(opp_hh));
					FacetHandle    fl = facet_handle(css);
					FacetHandle    fr = facet_handle(opp_hh);

					if (fl.is_valid()) {
						updateVertexPosition[i] += facet_ref(fl).normal_*(facet_ref(fl).normal_ *(opp_vc - vc));

					}
					if (fr.is_valid()) {
						updateVertexPosition[i] += facet_ref(fr).normal_*(facet_ref(fr).normal_ *(opp_vc - vc));

					}

					css = cw_rotated(css);

				} while (css != hh);
			}
			for (i = 0; i < vertex_num; i++) { vertex_ref(VertexHandle(i)).coord_ += updateVertexPosition[i] * 1.0 / 18.0; }

		} while (--iterations);

		clock_t t2 = clock();
		std::cout << "The time of Sun's vertex updating: " << (t2 - t1)*1.0 / CLOCKS_PER_SEC << "s" << std::endl;
	}

////////////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	void ExKernelT<ExItems>::mesh_process(){///////(��ѡ��������������ʵ��һ�ֲ���������������ȡ������ָ����������εȵ�)

		std::cout << "������˹��˳��������1����Ӹ�˹����������2\n";
		int i;
		std::cin >> i;
		switch (i) 
		{
		case 1:  //������˹��˳
		{
			int iterations;
			int i = 0;
			int vertex_num = vertex_size();
			std::cout << "�������������(10-30)��";
			std::cin >> iterations;
			std::cout << "���ڵ������ȽϺ�ʱ��" << std::endl;
			std::cout << "Please wait.....  " << std::endl;
			do {
				std::vector<Coord> updateVertexPosition;
				updateVertexPosition.resize(vertex_num);
				VertexIterator ff(vertex_begin());
				for (; ff != vertex_end(); ff++)
				{
					HalfedgeHandle& fh = ff->halfedge_handle_;
					VertexHandle& fv = vertex_handle(fh);
					HalfedgeHandle hhh(fh);                                          //�õ�����
					Coord coord(0, 0, 0);
					Scalar v = 0;                                               //�����������
					do {
						VertexHandle& fv = vertex_handle(opposite_halfedge_handle(hhh));//�ɰ����õ�������
						Coord& position = vertex_ref(fv).coord_;
						coord += position;
						v = v + 1.0;
						hhh = cw_rotated(hhh);
					} while (hhh != fh);
					coord = coord / v;
					updateVertexPosition[i++] = coord;
				}
				for (i = 0; i < vertex_num; i++)
					vertex_ref(VertexHandle(i)).coord_ = updateVertexPosition[i];//�������е�����
				i = 0;
			} while (iterations--);

		}
			break;
		case 2:   //��Ӹ�˹����
		{
			const double mean = 0.0;//��ֵ
			const double stddev = 0.01;//��׼��
			std::default_random_engine generator;   //��˹�������ɺ���
			std::normal_distribution<double> dist(mean, stddev);

			int vertex_num = vertex_size();
			std::vector<Coord> updateVertexPosition;
			updateVertexPosition.resize(vertex_num);
			for (i = 0; i < vertex_num; i++) {
				VertexHandle vh(i);
				Coord&       vc = coord(vh); //���ÿ�����������

				updateVertexPosition[i].data_[0] = vc.data_[0] + dist(generator);  //x������Ӹ�˹����
				updateVertexPosition[i].data_[1] = vc.data_[1] + dist(generator); //y������Ӹ�˹����
				updateVertexPosition[i].data_[2] = vc.data_[2] + dist(generator); //z������Ӹ�˹����
			}
			for (i = 0; i < vertex_num; i++) { vertex_ref(VertexHandle(i)).coord_ = updateVertexPosition[i]; }

		}
			break;
		default:
			std :: cout << "��������������������\n";
			break;
		}



	}

////////////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	typename ExKernelT<ExItems>::Scalar
		ExKernelT<ExItems>::calc_facet_area(const FacetHandle& _fh){/////���������ε����
		HalfedgeHandle& m = halfedge_handle(_fh);
		HalfedgeHandle& pre_m = prev_halfedge_handle(m);
		HalfedgeHandle& nex_m = next_halfedge_handle(m);    ///�õ���߾��

		EdgeHandle& n = edge_handle(m);
		EdgeHandle& pre_n = edge_handle(pre_m);
		EdgeHandle& nex_n = edge_handle(nex_m);      //�ɰ�߾���õ��߾��


		float a = calc_edge_length(n);
		float b = calc_edge_length(pre_n);
		float c = calc_edge_length(nex_n);         //�ɱ߾���õ����߳�
		float p = (a + b + c) / 2;

		Scalar 	area = sqrt(p*(p - a)*(p - b)*(p - c));////���ú��׹�ʽ�����s=sqrt(p(p-a)(p-b)(p-c))



		facet_ref(_fh).area_ = area;
		return area;
	}

	//////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	void ExKernelT<ExItems>::update_area()//��������������Ƭ�����
	{   
		set_isArea(true);
		Scalar max_area = 0.0;

		FacetIterator fit(facet_begin() );
		for(;fit<facet_end(); fit++)
		{
			HalfedgeHandle hh = fit->halfedge_handle_;
			FacetHandle    fh = facet_handle(hh);
			Scalar curr_area = calc_facet_area(fh);

			if (curr_area > max_area)
				max_area = curr_area;
		}

		std::cout << "The maximal area of the mesh is: " << max_area << std::endl;
	}
	///////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	typename ExKernelT<ExItems>::Scalar
		ExKernelT<ExItems>::get_area(const FacetHandle& _fh)
	{
		return facet_ref(_fh).area_;
	}
	//////////////////////////////////////////////////////////////////////////////
	template<class ExItems> 
	double 
		ExKernelT<ExItems>::calc_edge_length(const EdgeHandle &_eh) {//updating an edge length


			HalfedgeHandle& h1 = halfedge_handle(_eh,0);
			HalfedgeHandle& h2 = halfedge_handle(_eh,1);

			Vertex v0 = vertex_ref(vertex_handle(h1) );
			Vertex v1 = vertex_ref(vertex_handle(h2) );

			return (v0.coord_-v1.coord_).norm();

	}

	template<class ExItems>
	void
		ExKernelT<ExItems>::getNeighborRing(VertexHandle& _vh, int _ring, std::vector<VertexHandle>& NeighborRing){//�õ���ring���ڽӵ�

			NeighborRing.push_back( _vh );
			int iteration = 0;
			int verNewNum = NeighborRing.size();
			int verOldNum = verNewNum-1;
			int verOldNum1 = verOldNum;

			do{
				verOldNum = NeighborRing.size();
				for(int i=verOldNum1; i<verNewNum; i++){
					VertexHandle& vh = NeighborRing[i];
					HalfedgeHandle& hh = halfedge_handle(vh);
					HalfedgeHandle css(hh);
					do{
						HalfedgeHandle& opp_hh = opposite_halfedge_handle(css);
						VertexHandle&   opp_vh = vertex_handle(opp_hh);
						for(int ii=0; ii<NeighborRing.size(); ii++)
							if(opp_vh == NeighborRing[ii] ) break;

						if(ii >= NeighborRing.size() )
							NeighborRing.push_back(opp_vh);

						css = cw_rotated(css);
					}while(css != hh);
				}

				verNewNum = NeighborRing.size();
				verOldNum1 = verOldNum;
				iteration++;
			}while(iteration < _ring);

	}
	///////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	void
		ExKernelT<ExItems>::getNeighborFaceN1(FacetHandle& _fh, std::vector<FacetHandle>& _fhs){//sharing common edges

			HalfedgeHandle& hh = halfedge_handle(_fh);
			HalfedgeHandle& pre_hh = prev_halfedge_handle(hh);
			HalfedgeHandle& nex_hh = next_halfedge_handle(hh);

			_fhs.push_back( facet_handle( opposite_halfedge_handle(hh) ) );
			_fhs.push_back( facet_handle( opposite_halfedge_handle(pre_hh ) ) );
			_fhs.push_back( facet_handle(opposite_halfedge_handle(nex_hh ) ) );
	}
	///////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	void 
		ExKernelT<ExItems>::getNeighborFaceN2(FacetHandle& _fh, std::vector<FacetHandle>& _fhs){//sharing common vertices

			HalfedgeHandle& hh = halfedge_handle(_fh);
			HalfedgeHandle& pre_hh = prev_halfedge_handle(hh);
			HalfedgeHandle& nex_hh = next_halfedge_handle(hh);

			VertexHandle  vhs[3];
			vhs[0] = vertex_handle(hh);
			vhs[1] = vertex_handle(pre_hh);
			vhs[2] = vertex_handle(nex_hh);
			int i = 0, j=0;

			for(i=0; i<3; i++){

				HalfedgeHandle& hhv = halfedge_handle( vhs[i] );
				HalfedgeHandle cursor(hhv);

				do{

					FacetHandle fh = facet_handle(cursor);
					if(fh.is_valid() && fh != _fh){

						if(_fhs.size() != 0){

							for(j=0; j< _fhs.size(); j++){

								if(fh.idx() == _fhs[j].idx() ) break;
							}//end for

							if(j>= _fhs.size() ) _fhs.push_back(fh);

						}//end if
						else _fhs.push_back(fh);
					}//end if

					cursor = cw_rotated(cursor);

				}while(hhv != cursor);//end for do while
			}//end for

	}
	///////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	void ExKernelT<ExItems>::output_to_file(){

		FILE *fp;
		fp=fopen("My_result.off","w");

		int no_vertex=vertex_size();
		int no_facet=facet_size();
		int edge = 0;

		fprintf(fp,"OFF\n");
		fprintf(fp,"%d  %d %d\n",no_vertex,no_facet, edge);

		VertexIterator vit(vertex_begin());

		for(;vit!=vertex_end();vit++){

			VertexHandle vh=vertex_handle(vit->halfedge_handle_);
			fprintf(fp," %f  %f  %f\n",coord(vh).data_[0],coord(vh).data_[1],coord(vh).data_[2]);

		}

		FacetIterator fit(facet_begin());

		for(;fit!=facet_end();fit++){

			HalfedgeHandle hh=fit->halfedge_handle_;
			HalfedgeHandle nh=next_halfedge_handle(hh);
			HalfedgeHandle nnh=next_halfedge_handle(nh);

			VertexHandle vh=vertex_handle(hh);
			VertexHandle nvh=vertex_handle(nh);
			VertexHandle nnvh=vertex_handle(nnh);

			fprintf(fp,"3 %d  %d  %d\n",vh.idx(),nvh.idx(),nnvh.idx());

		}

		fclose(fp);
	}
//////////////////////////////////////////////////////////////////////////////////
	template <class ExItems> 
	typename ExKernelT<ExItems>::Normal
		ExKernelT<ExItems>::normal(const FacetHandle& _fh) {
			assert( _fh.is_valid() );
			assert( _fh.idx() < facet_size() );
			return facet_ref(_fh).normal_;
	}


	///////////////////////////////////////////////////////////////////////////////
	template <class ExItems> 
	typename ExKernelT<ExItems>::Normal
		ExKernelT<ExItems>::calc_normal(const FacetHandle& _fh) {
			assert( _fh.is_valid() );
			assert( _fh.idx() < facet_size() );

			const HalfedgeHandle&   hh = halfedge_handle(_fh);
			const HalfedgeHandle& p_hh = prev_halfedge_handle(hh);
			const HalfedgeHandle& n_hh = next_halfedge_handle(hh);

			const Coord& cd0 = coord( vertex_handle( hh) );
			const Coord& cd1 = coord( vertex_handle(p_hh) );
			const Coord& cd2 = coord( vertex_handle(n_hh) );

			//return ((cd1-cd0)%(cd2-cd1)).normalize();
			return ((cd2-cd1)%(cd1-cd0)).normalize();//����������˲��ҵ�λ�� be careful
	}


	///////////////////////////////////////////////////////////////////////////////
	template <class ExItems> 
	void ExKernelT<ExItems>::update_facet_normals(void) {

		set_isNormal(true);
		FacetIterator fi = facet_begin();

		for ( ; fi!=facet_end(); ++fi) {
			if (fi->status_.is_deleted()) continue;

			assert(fi->halfedge_handle_.is_valid());
			fi->normal_ = calc_normal( facet_handle(fi->halfedge_handle_) );
		}
	}


	///////////////////////////////////////////////////////////////////////////////
	template <class ExItems> 
	typename ExKernelT<ExItems>::Normal
		ExKernelT<ExItems>::normal(const VertexHandle& _vh) {
			assert( _vh.is_valid() );
			assert( _vh.idx() < vertex_size() );
			return vertex_ref(_vh).normal_;
	}


///////////////////////////////////////////////////////////////////////////////
	template <class ExItems> 
	typename ExKernelT<ExItems>::Normal
		ExKernelT<ExItems>::calc_normal(const VertexHandle& _vh) {////����һ������ķ���
		assert(_vh.is_valid());
		assert(_vh.idx() < vertex_size());

		Normal          norm(0, 0, 0);///////��norm�洢��õķ���ֵ

//////////////////////////�ڴ�ʵ��////////////////////////////////////
		double allarea = 0;
		HalfedgeHandle& hhh = halfedge_handle(_vh);
		HalfedgeHandle csss(hhh);
		do {                                           //�����������������Ƭ�����������Ƭ�������allarea
			FacetHandle& dd = facet_handle(csss);
			if (dd.is_valid())
				allarea += calc_facet_area(dd);
			csss = cw_rotated(csss);
		} while (csss != hhh);

		HalfedgeHandle& hh = halfedge_handle(_vh);
		HalfedgeHandle css(hh);
		do {
			FacetHandle& d = facet_handle(css);
			if (d.is_valid())
				norm += calc_normal(d)*(calc_facet_area(d) / allarea);           //������Ƭ�����Ȩ������Ƭ���� * ����������Ƭ�����Ȩֵ
			css = cw_rotated(css);
		} while (css != hh);


		/////////////////////////////////////////////////////////////////////

		return norm.normalize();  //��λ��
	}


///////////////////////////////////////////////////////////////////////////////
	template <class ExItems> 
	void ExKernelT<ExItems>::update_vertex_normals(void) {//�������ж���ķ���
		VertexIterator vi = vertex_begin();

		for ( ; vi!=vertex_end(); ++vi) {
			if (vi->status_.is_deleted()) continue;

			assert(vi->halfedge_handle_.is_valid());
			vi->normal_ = calc_normal( vertex_handle(vi->halfedge_handle_) );
		}
	}

///////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	typename ExKernelT<ExItems>::Normal
		ExKernelT<ExItems>::calc_normal_max(const VertexHandle& _vh){

			assert(_vh.is_valid() );
			assert(_vh.idx() < vertex_size() );

			HalfedgeHandle& hh = halfedge_handle(_vh);
			Coord& vc = vertex_ref(_vh).coord_;
			Normal n(0,0,0);
			HalfedgeHandle css(hh);
			do{
				FacetHandle& fh = facet_handle(css);
				if(fh.is_valid() ){

					HalfedgeHandle& nexhh = next_halfedge_handle(css);
					HalfedgeHandle& prehh = prev_halfedge_handle(css);

					VertexHandle& nvh = vertex_handle(nexhh);
					VertexHandle& prevh = vertex_handle(prehh);

					Coord& nvc = vertex_ref(nvh).coord_;
					Coord& prevc = vertex_ref(prevh).coord_;

					Coord vec1 = vc - nvc;
					Coord vec2 = vc - prevc;
					Coord vec12cross = vec1 % vec2;//cross multiplication

					float weight = vec12cross.length() / (vec1.sqLength() * vec2.sqLength() );

					n += facet_ref(fh).normal_ * weight;//---------------------------------------ע��
					//n += calc_normal(fh);//---------------------------------ע��

				}//if

				css = cw_rotated(css);
			}while(css != hh);

			n.normalize();

			return n;

	}
	//////////////////////////////////////////////////////////////////////////////////
	template <class ExItems> 
	void ExKernelT<ExItems>::update_vertex_normals_max(void){

		VertexIterator vi = vertex_begin();

		for ( ; vi!=vertex_end(); ++vi) {
			if (vi->status_.is_deleted() ) continue;

			assert(vi->halfedge_handle_.is_valid());
			//vi->normal_ = calc_normal( vertex_handle(vi->halfedge_handle_));
			vi->normal_ = calc_normal_max( vertex_handle(vi->halfedge_handle_) );
			//std::cout<<vi->normal_<<std::endl;
		}
	}
	///////////////////////////////////////////////////////////////////////////////////
	template <class ExItems> 
	void ExKernelT<ExItems>::update_normals(void) {
		update_facet_normals();
		update_vertex_normals();
	}


	///////////////////////////////////////////////////////////////////////////////
	template <class ExItems> 
	void ExKernelT<ExItems>::update_edge_length(void) {//���������б߳���Ϣ
		float global_max_edge_length_ = 0;
		float averagedlength = 0.0;

		EdgeIterator eit = edge_begin();
		for ( ; eit!=edge_end(); ++eit ) {
			Vertex& v0 = vertex_ref( eit->halfedges_[0].vertex_handle_ );
			Vertex& v1 = vertex_ref( eit->halfedges_[1].vertex_handle_ );

			eit->length_ = (v0.coord_ - v1.coord_).norm();
			//std::cout<<eit->length_<<"\n";
			averagedlength += eit->length_;

			if (global_max_edge_length_ < eit->length_)
				global_max_edge_length_ = eit->length_; 
		} 

		averagedlength /= edge_size();
		set_average_edge_length(averagedlength);
	}
	/////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	typename ExKernelT<ExItems>::Coord
		ExKernelT<ExItems>::calc_centroid(const FacetHandle& _fh){

			HalfedgeHandle& hh = halfedge_handle(_fh);
			HalfedgeHandle& n_hh = next_halfedge_handle(hh);
			HalfedgeHandle& pre_hh = prev_halfedge_handle(hh);

			VertexHandle& vh = vertex_handle(hh);
			VertexHandle& n_vh = vertex_handle(n_hh);
			VertexHandle& pre_vh = vertex_handle(pre_hh);

			return Coord(coord(vh) +
				coord(n_vh)+
				coord(pre_vh) )/3.0;

	}

	//////////////////////////////////////////////////////////////////////////
	//template<class ExItems>
	//bool ExKernelT<ExItems>::createPS(){

	//	unsigned int NumFacet = facet_size();

	//	ps_ = new PointSet(NumFacet);//��ʼ��

	//	FacetIterator fi = facet_begin();
	//	unsigned int i = 0;

	//	for(; fi != facet_end(); fi++ ){

	//		Coord& cor = calc_centroid( facet_handle(fi->halfedge_handle_)
	//			);
	//		ps_->setPoint(i, cor[0], cor[1], cor[2] );

	//		i++;
	//	}

	//	//std::cout<<ps_->pointN<<"pointN"<<std::endl;

	//	return true;

	//}
	///////////////////////////////////////////////////////////////
	//template<class ExItems>
	//bool ExKernelT<ExItems>::createKdTree(){

	//	kdtree_ = new KdTree(ps_);

	//	//std::cout<<kdtree_->listN<<"listN"<<std::endl;

	//	return true;


	//}
	
////////////////////////////////////////////////////////////////
	template<class ExItems>
	void
		ExKernelT<ExItems>::adjustFaceNormal_FE(){///Sun Xiangfang TVCG 2007 fast and effective feature-preserving mesh denoising

			int facet_num = facet_size();

			double T;
			std::cout<<"Input Threshold(0.0-1.0): ";
			std::cin>>T;

			int iteration;
			std::cout<<"Input facenormal filtering numbers (5-20��): ";
			std::cin>>iteration;
			std::cout << "���ڵ������ȽϺ�ʱ��"<<std::endl;
			std::cout << "Please wait.....  "<<std::endl;
			int i = 0, j = 0;
			clock_t t1 = clock();
			do{

				std::vector<Normal> updateFacetNormal;
				updateFacetNormal.resize(facet_num);

				FacetIterator fit = facet_begin();
				for( ; fit!=facet_end(); fit++){
					std::vector<FacetHandle> fhs;
					HalfedgeHandle hh = fit->halfedge_handle_;
					FacetHandle    fh = facet_handle(hh);
					fhs.push_back(fh);
					//VertexHandle varry[3];
					std::vector<VertexHandle> varry;
					varry.resize(3);
					varry[0] = vertex_handle(hh);
					varry[1] = vertex_handle( prev_halfedge_handle(hh) );
					varry[2] = vertex_handle( next_halfedge_handle(hh) );
					for(i=0; i<3; i++){
						HalfedgeHandle vahh = halfedge_handle(varry[i] );
						HalfedgeHandle vacss(vahh);
						do{
							FacetHandle vafh = facet_handle(vacss );
							if(vafh.is_valid() ){

								for(j=0;j<fhs.size(); j++)
									if(vafh == fhs[j] ) break;
								if(j>=fhs.size() ) fhs.push_back(vafh);
							}

							vacss = cw_rotated(vacss );

						}while(vacss != vahh);
					}

					Normal& nf = facet_ref(fh).normal_;
					for(i=0; i<fhs.size(); i++){

						Normal& nfi = facet_ref(fhs[i]).normal_;
						double dotproduct = nf*nfi - T;
						//double dotproduct = 1 - nf*nfi;
						if(dotproduct <1e-7) continue;
						//dotproduct = exp(dotproduct);

						dotproduct *= dotproduct;

						updateFacetNormal[fh.idx()] += nfi * dotproduct;


					}


				}

				for(fit=facet_begin();fit != facet_end(); fit++){
					facet_ref(facet_handle(fit->halfedge_handle_ ) ).normal_ = updateFacetNormal[facet_handle(fit->halfedge_handle_).idx()].normalize();
				}

			}while(--iteration);

			clock_t t2 = clock();

			//t2 = (t2-t1)*1.0/CLOCKS_PER_SEC;
			std::cout<<"The time of Sun's normal filter: "<<(t2-t1)*1.0/CLOCKS_PER_SEC<<"s"<<std::endl;

	}
	////////////////////////////////////////////////////////////////
	template<class ExItems>
	void 
		ExKernelT<ExItems>::Mesh_Denoising_FE(){//Sun Xianfang TVCG2007
			adjustFaceNormal_FE();
			//adjustFaceNormal_Based_on_FE();

			int vertex_num = vertex_size();
			int iterations;
			std::cout<<"Input vertex update iterations (10-30��): ";
			std::cin>>iterations;
			std::cout << "���ڵ������ȽϺ�ʱ��" << std::endl;
			std::cout << "Please wait.....  " << std::endl;
			int i = 0;

			clock_t t1 = clock();
			do{
				std::vector<Coord> updateVertexPosition;
				updateVertexPosition.resize(vertex_num);
				for(i=0; i<vertex_num; i++){
					VertexHandle vh(i);
					Coord&       vc = coord(vh);
					HalfedgeHandle& hh = halfedge_handle(vh);
					HalfedgeHandle  css(hh);
					do{
						HalfedgeHandle opp_hh = opposite_halfedge_handle(css);
						Coord&         opp_vc = coord(vertex_handle(opp_hh) );
						FacetHandle    fl = facet_handle(css);
						FacetHandle    fr = facet_handle(opp_hh);

						if(fl.is_valid() ){
							updateVertexPosition[i] += facet_ref(fl).normal_*( facet_ref(fl).normal_ *(opp_vc - vc) );

						}
						if(fr.is_valid() ){
							updateVertexPosition[i] += facet_ref(fr).normal_*( facet_ref(fr).normal_ *(opp_vc - vc) );

						}

						css = cw_rotated(css);

					}while(css != hh);
				}
				for(i=0; i<vertex_num; i++){ vertex_ref(VertexHandle(i) ).coord_ += updateVertexPosition[i]*1.0/18.0;}

			}while(--iterations);

			clock_t t2 = clock();
			//t2 = (t2-t1)/CLOCKS_PER_SEC;
			std::cout<<"The time of Sun's vertex updating: "<<(t2-t1)*1.0/CLOCKS_PER_SEC<<"s"<<std::endl;

	}
	////////////////////////////////////////////////
	template<class ExItems>
	void 
		ExKernelT<ExItems>::adjustFaceNormal_YouyiZheng(){

			double SigmaC = get_average_edge_length();
			float MinSigmaS, MaxSigmaS;
			//calcdSigmaS1(MinSigmaS,MaxSigmaS);
			//MaxSigmaS /= 2.0;
			std::cout<<"Please input SigmaS (0.2-0.6): ";
			std::cin>>MaxSigmaS;

			std::vector<double> facetareas;//getting each facet area
			FacetIterator fi(facet_begin() );
			for(; fi < facet_end(); fi++){

				HalfedgeHandle& hh = fi->halfedge_handle_;
				FacetHandle&    fh = facet_handle(hh);
				facetareas.push_back(calc_facet_area(fh) );

			}//for

			std::cout<<"Input facet normal filtering iteration numbers (5-20��): ";
			int num;//the iteration num
			std::cin>>num;
			std::cout << "���ڵ������ȽϺ�ʱ��" << std::endl;
			std::cout << "Please wait.....  " << std::endl;

			clock_t t1 = clock();
			int i = 0;

			for(i=0; i<num; i++){

				std::vector<Normal> new_facet_normals;

				for(fi = facet_begin(); fi < facet_end(); fi++){
					HalfedgeHandle& hh = fi->halfedge_handle_;
					FacetHandle&    fh = facet_handle(hh);
					FacetHandles    NeighorFacets;
					//getNeighborFaceN1(fh,NeighorFacets); //getting facets sharing common edges
					getNeighborFaceN2(fh,NeighorFacets);//sharing common vertices

					Coord& center_fh = calc_centroid(fh);
					Normal& nf = normal(fh);
					Normal  trans(0,0,0);
					for(int j=0; j<NeighorFacets.size(); j++){

						FacetHandle& Nfh = NeighorFacets[j];
						Coord& Ncenter = calc_centroid(Nfh);
						Normal& NorN   = normal(Nfh);

						double Wc, Ws;
						Wc = exp( -(center_fh - Ncenter).sqNorm()/ (2*SigmaC*SigmaC) );
						Ws = exp( -(nf - NorN).sqNorm()          / (2*MaxSigmaS*MaxSigmaS) );

						trans += NorN * Wc*Ws*facetareas[Nfh.idx() ];


					}// end for j

					if( trans.length() > 1E-7) new_facet_normals.push_back( trans.normalize() );
					else new_facet_normals.push_back(trans);

				}//for end fi

				for(fi=facet_begin(); fi<facet_end(); fi++){

					FacetHandle fhh = facet_handle(fi->halfedge_handle_);
					fi->normal_ = new_facet_normals[fhh.idx() ];
				}
			}//for end i

			clock_t t2 = clock();
			//t2 = (t2-t1)/CLOCKS_PER_SEC;
			std::cout<<"The time of Zheng's normal filter: "<<(t2-t1)*1.0/CLOCKS_PER_SEC<<"s"<<std::endl;

	}
	///////////////////////////////////////////////////////////////
	template<class ExItems>
	void 
		ExKernelT<ExItems>::Mesh_Denoising_YouyiZheng(){//according to YouyiZheng TVCG 2011

			adjustFaceNormal_YouyiZheng();

			int vertex_num = vertex_size();
			int iterations;
			std::cout<<"Input vertex update iterations(5-30��): ";
			std::cin>>iterations;
			int i = 0;
			std::cout << "���ڵ������ȽϺ�ʱ��" << std::endl;
			std::cout << "Please wait.....  " << std::endl;
			clock_t t1 = clock();
			do{
				std::vector<Coord> updateVertexPosition;
				updateVertexPosition.resize(vertex_num);
				for(i=0; i<vertex_num; i++){
					VertexHandle vh(i);
					Coord&       vc = coord(vh);
					HalfedgeHandle& hh = halfedge_handle(vh);
					HalfedgeHandle  css(hh);
					do{
						HalfedgeHandle opp_hh = opposite_halfedge_handle(css);
						Coord&         opp_vc = coord(vertex_handle(opp_hh) );
						FacetHandle    fl = facet_handle(css);
						FacetHandle    fr = facet_handle(opp_hh);

						if(fl.is_valid() ){
							updateVertexPosition[i] += facet_ref(fl).normal_*( facet_ref(fl).normal_ *(opp_vc - vc) );

						}
						if(fr.is_valid() ){
							updateVertexPosition[i] += facet_ref(fr).normal_*( facet_ref(fr).normal_ *(opp_vc - vc) );

						}

						css = cw_rotated(css);

					}while(css != hh);
				}
				for(i=0; i<vertex_num; i++){ vertex_ref(VertexHandle(i) ).coord_ += updateVertexPosition[i]*1.0/18.0;}

			}while(--iterations);

			clock_t t2 = clock();
			//t2 = (t2-t1)/CLOCKS_PER_SEC;
			std::cout<<"The time of Zheng's vertex updating: "<<(t2-t1)*1.0/CLOCKS_PER_SEC<<"s"<<std::endl;

	}
	///////////////////////////////////////////////////////////////

} /// namespace




