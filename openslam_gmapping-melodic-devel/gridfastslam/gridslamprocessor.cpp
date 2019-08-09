#include <string>
#include <deque>
#include <list>
#include <map>
#include <set>
#include <fstream>
#include <iomanip>
#include <gmapping/utils/stat.h>
#include "gmapping/gridfastslam/gridslamprocessor.h"

//#define MAP_CONSISTENCY_CHECK
//#define GENERATE_TRAJECTORIES

namespace GMapping {

const double m_distanceThresholdCheck = 20;
 
using namespace std;

  GridSlamProcessor::GridSlamProcessor(): m_infoStream(cout){
    
    period_ = 5.0;
    m_obsSigmaGain=1;
    m_resampleThreshold=0.5;
    m_minimumScore=0.;
  }
  
  GridSlamProcessor::GridSlamProcessor(const GridSlamProcessor& gsp) 
    :last_update_time_(0.0), m_particles(gsp.m_particles), m_infoStream(cout){
    
    period_ = 5.0;

    m_obsSigmaGain=gsp.m_obsSigmaGain;
    m_resampleThreshold=gsp.m_resampleThreshold;
    m_minimumScore=gsp.m_minimumScore;
    
    m_beams=gsp.m_beams;
    m_indexes=gsp.m_indexes;
    m_motionModel=gsp.m_motionModel;
    m_resampleThreshold=gsp.m_resampleThreshold;
    m_matcher=gsp.m_matcher;
    
    m_count=gsp.m_count;
    m_readingCount=gsp.m_readingCount;
    m_lastPartPose=gsp.m_lastPartPose;
    m_pose=gsp.m_pose;
    m_odoPose=gsp.m_odoPose;
    m_linearDistance=gsp.m_linearDistance;
    m_angularDistance=gsp.m_angularDistance;
    m_neff=gsp.m_neff;
	
    cerr << "FILTER COPY CONSTRUCTOR" << endl;
    cerr << "m_odoPose=" << m_odoPose.x << " " <<m_odoPose.y << " " << m_odoPose.theta << endl;
    cerr << "m_lastPartPose=" << m_lastPartPose.x << " " <<m_lastPartPose.y << " " << m_lastPartPose.theta << endl;
    cerr << "m_linearDistance=" << m_linearDistance << endl;
    cerr << "m_angularDistance=" << m_linearDistance << endl;
    
		
    m_xmin=gsp.m_xmin;
    m_ymin=gsp.m_ymin;
    m_xmax=gsp.m_xmax;
    m_ymax=gsp.m_ymax;
    m_delta=gsp.m_delta;
    
    m_regScore=gsp.m_regScore;
    m_critScore=gsp.m_critScore;
    m_maxMove=gsp.m_maxMove;
    
    m_linearThresholdDistance=gsp.m_linearThresholdDistance;
    m_angularThresholdDistance=gsp.m_angularThresholdDistance;
    m_obsSigmaGain=gsp.m_obsSigmaGain;
    
#ifdef MAP_CONSISTENCY_CHECK
    cerr << __PRETTY_FUNCTION__ <<  ": trajectories copy.... ";
#endif
    TNodeVector v=gsp.getTrajectories();
    for (unsigned int i=0; i<v.size(); i++){
		m_particles[i].node=v[i];
    }
#ifdef MAP_CONSISTENCY_CHECK
    cerr <<  "end" << endl;
#endif


    cerr  << "Tree: normalizing, resetting and propagating weights within copy construction/cloneing ..." ;
    updateTreeWeights(false);
    cerr  << ".done!" <<endl;
  }
  
  GridSlamProcessor::GridSlamProcessor(std::ostream& infoS): m_infoStream(infoS){
    period_ = 5.0;
    m_obsSigmaGain=1;
    m_resampleThreshold=0.5;
    m_minimumScore=0.;
	
  }

  GridSlamProcessor* GridSlamProcessor::clone() const {
# ifdef MAP_CONSISTENCY_CHECK
    cerr << __PRETTY_FUNCTION__ << ": performing preclone_fit_test" << endl;
    typedef std::map<autoptr< Array2D<PointAccumulator> >::reference* const, int> PointerMap;
    PointerMap pmap;
	for (ParticleVector::const_iterator it=m_particles.begin(); it!=m_particles.end(); it++){
	  const ScanMatcherMap& m1(it->map);
	  const HierarchicalArray2D<PointAccumulator>& h1(m1.storage());
 	  for (int x=0; x<h1.getXSize(); x++){
	    for (int y=0; y<h1.getYSize(); y++){
	      const autoptr< Array2D<PointAccumulator> >& a1(h1.m_cells[x][y]);
	      if (a1.m_reference){
		PointerMap::iterator f=pmap.find(a1.m_reference);
		if (f==pmap.end())
		  pmap.insert(make_pair(a1.m_reference, 1));
		else
		  f->second++;
	      }
	    }
	  }
	}
	cerr << __PRETTY_FUNCTION__ <<  ": Number of allocated chunks" << pmap.size() << endl;
	for(PointerMap::const_iterator it=pmap.begin(); it!=pmap.end(); it++)
	  assert(it->first->shares==(unsigned int)it->second);

	cerr << __PRETTY_FUNCTION__ <<  ": SUCCESS, the error is somewhere else" << endl;
# endif
	GridSlamProcessor* cloned=new GridSlamProcessor(*this);
	
# ifdef MAP_CONSISTENCY_CHECK
	cerr << __PRETTY_FUNCTION__ <<  ": trajectories end" << endl;
	cerr << __PRETTY_FUNCTION__ << ": performing afterclone_fit_test" << endl;
	ParticleVector::const_iterator jt=cloned->m_particles.begin();
	for (ParticleVector::const_iterator it=m_particles.begin(); it!=m_particles.end(); it++){
	  const ScanMatcherMap& m1(it->map);
	  const ScanMatcherMap& m2(jt->map);
	  const HierarchicalArray2D<PointAccumulator>& h1(m1.storage());
	  const HierarchicalArray2D<PointAccumulator>& h2(m2.storage());
	  jt++;
 	  for (int x=0; x<h1.getXSize(); x++){
	    for (int y=0; y<h1.getYSize(); y++){
	      const autoptr< Array2D<PointAccumulator> >& a1(h1.m_cells[x][y]);
	      const autoptr< Array2D<PointAccumulator> >& a2(h2.m_cells[x][y]);
	      assert(a1.m_reference==a2.m_reference);
	      assert((!a1.m_reference) || !(a1.m_reference->shares%2));
	    }
	  }
	}
	cerr << __PRETTY_FUNCTION__ <<  ": SUCCESS, the error is somewhere else" << endl;
# endif
	return cloned;
}
  
  GridSlamProcessor::~GridSlamProcessor(){
    cerr << __PRETTY_FUNCTION__ << ": Start" << endl;
    cerr << __PRETTY_FUNCTION__ << ": Deleting tree" << endl;
    for (std::vector<Particle>::iterator it=m_particles.begin(); it!=m_particles.end(); it++){
#ifdef TREE_CONSISTENCY_CHECK		
      TNode* node=it->node;
      while(node)
	node=node->parent;
      cerr << "@" << endl;
#endif
      if (it->node)
	delete it->node;
      //cout << "l=" << it->weight<< endl;
    }
    
# ifdef MAP_CONSISTENCY_CHECK
    cerr << __PRETTY_FUNCTION__ << ": performing predestruction_fit_test" << endl;
    typedef std::map<autoptr< Array2D<PointAccumulator> >::reference* const, int> PointerMap;
    PointerMap pmap;
    for (ParticleVector::const_iterator it=m_particles.begin(); it!=m_particles.end(); it++){
      const ScanMatcherMap& m1(it->map);
      const HierarchicalArray2D<PointAccumulator>& h1(m1.storage());
      for (int x=0; x<h1.getXSize(); x++){
	for (int y=0; y<h1.getYSize(); y++){
	  const autoptr< Array2D<PointAccumulator> >& a1(h1.m_cells[x][y]);
	  if (a1.m_reference){
	    PointerMap::iterator f=pmap.find(a1.m_reference);
	    if (f==pmap.end())
	      pmap.insert(make_pair(a1.m_reference, 1));
	    else
	      f->second++;
	  }
	}
      }
    }
    cerr << __PRETTY_FUNCTION__ << ": Number of allocated chunks" << pmap.size() << endl;
    for(PointerMap::const_iterator it=pmap.begin(); it!=pmap.end(); it++)
      assert(it->first->shares>=(unsigned int)it->second);
    cerr << __PRETTY_FUNCTION__ << ": SUCCESS, the error is somewhere else" << endl;
# endif
  }


		
  void GridSlamProcessor::setMatchingParameters (double urange, double range, double sigma, int kernsize, double lopt, double aopt, 
						 int iterations, double likelihoodSigma, double likelihoodGain, unsigned int likelihoodSkip){
    m_obsSigmaGain=likelihoodGain;
    m_matcher.setMatchingParameters(urange, range, sigma, kernsize, lopt, aopt, iterations, likelihoodSigma, likelihoodSkip);
    if (m_infoStream)
      m_infoStream << " -maxUrange "<< urange
		   << " -maxUrange "<< range
		   << " -sigma     "<< sigma
		   << " -kernelSize "<< kernsize
		   << " -lstep "    << lopt
		   << " -lobsGain " << m_obsSigmaGain
		   << " -astep "    << aopt << endl;
    
    
  }
  
void GridSlamProcessor::setMotionModelParameters
(double srr, double srt, double str, double stt){
  m_motionModel.srr=srr;
  m_motionModel.srt=srt;
  m_motionModel.str=str;
  m_motionModel.stt=stt;	
  
  if (m_infoStream)
    m_infoStream << " -srr "<< srr 	<< " -srt "<< srt  
		 << " -str "<< str 	<< " -stt "<< stt << endl;
  
}
  
  void GridSlamProcessor::setUpdateDistances(double linear, double angular, double resampleThreshold){
    m_linearThresholdDistance=linear; 
    m_angularThresholdDistance=angular;
    m_resampleThreshold=resampleThreshold;	
    if (m_infoStream)
      m_infoStream << " -linearUpdate " << linear
		   << " -angularUpdate "<< angular
		   << " -resampleThreshold " << m_resampleThreshold << endl;
  }
  
  //HERE STARTS THE BEEF

  GridSlamProcessor::Particle::Particle(const ScanMatcherMap& m):
    map(m), pose(0,0,0), weight(0), weightSum(0), gweight(0), previousIndex(0){
    node=0;
  }
  
  
  void GridSlamProcessor::setSensorMap(const SensorMap& smap){
    
    /*
      Construct the angle table for the sensor
      
      FIXME For now detect the readings of only the front laser, and assume its pose is in the center of the robot 
    */
    
    SensorMap::const_iterator laser_it=smap.find(std::string("FLASER"));
    if (laser_it==smap.end()){
      cerr << "Attempting to load the new carmen log format" << endl;
      laser_it=smap.find(std::string("ROBOTLASER1"));
      assert(laser_it!=smap.end());
    }
    const RangeSensor* rangeSensor=dynamic_cast<const RangeSensor*>((laser_it->second));
    assert(rangeSensor && rangeSensor->beams().size());
    
    m_beams=static_cast<unsigned int>(rangeSensor->beams().size());
    double* angles=new double[rangeSensor->beams().size()];
    for (unsigned int i=0; i<m_beams; i++){
      angles[i]=rangeSensor->beams()[i].pose.theta;
    }
    m_matcher.setLaserParameters(m_beams, angles, rangeSensor->getPose());
    delete [] angles;
  }
  
  void GridSlamProcessor::init(unsigned int size, double xmin, double ymin, double xmax, double ymax, double delta, OrientedPoint initialPose){
    m_xmin=xmin;
    m_ymin=ymin;
    m_xmax=xmax;
    m_ymax=ymax;
    m_delta=delta;
    if (m_infoStream)
      m_infoStream 
	<< " -xmin "<< m_xmin
	<< " -xmax "<< m_xmax
	<< " -ymin "<< m_ymin
	<< " -ymax "<< m_ymax
	<< " -delta "<< m_delta
	<< " -particles "<< size << endl;
    

    m_particles.clear();
    TNode* node=new TNode(initialPose, 0, 0, 0);
    ScanMatcherMap lmap(Point(xmin+xmax, ymin+ymax)*.5, xmax-xmin, ymax-ymin, delta);
    for (unsigned int i=0; i<size; i++){
      m_particles.push_back(Particle(lmap));
      m_particles.back().pose=initialPose;
      m_particles.back().previousPose=initialPose;
      m_particles.back().setWeight(0);
      m_particles.back().previousIndex=0;
      
		// this is not needed
		//		m_particles.back().node=new TNode(initialPose, 0, node, 0);

		// we use the root directly
		m_particles.back().node= node;
    }
    m_neff=(double)size;
    m_count=0;
    m_readingCount=0;
    m_linearDistance=m_angularDistance=0;
  }

  void GridSlamProcessor::processTruePos(const OdometryReading& o){
    const OdometrySensor* os=dynamic_cast<const OdometrySensor*>(o.getSensor());
    if (os && os->isIdeal() && m_outputStream){
      m_outputStream << setiosflags(ios::fixed) << setprecision(3);
      m_outputStream <<  "SIMULATOR_POS " <<  o.getPose().x << " " << o.getPose().y << " " ;
      m_outputStream << setiosflags(ios::fixed) << setprecision(6) << o.getPose().theta << " " <<  o.getTime() << endl;
    }
  }

/* 
双线程和程序的基本执行流程

GMapping采用GridSlamProcessorThread后台线程进行数据的读取和运算，在gsp_thread.h和相应的实现文件
gsp_thread.cpp中可以看到初始化，参数配置，扫描数据读取等方法。

最关键的流程是在后台线程调用的方法

void * GridSlamProcessorThread::fastslamthread(GridSlamProcessorThread* gpt)

而这个方法中最主要的处理扫描数据的过程在428行，
bool processed=gpt->processScan(*rr);同时可以看到GMapping支持在线和离线两种模式。

注意到struct GridSlamProcessorThread : public GridSlamProcessor ，
这个后台线程执行类GridSlamProcessorThread 继承自GridSlamProcessor，
所以执行的是GridSlamProcessor的processScan方法。

*/
  
  //后台线程处理扫描数据方法
  bool GridSlamProcessor::processScan(const RangeReading & reading, int adaptParticles){
    
    /**retireve the position from the reading, and compute the odometry*/

	/* 得到当前的里程计的位置        */
	
    OrientedPoint relPose=reading.getPose();

	/*m_count表示这个函数被调用的次数 如果是第0次调用,则所有的位姿都是一样的*/
	
    if (!m_count){
      m_lastPartPose=m_odoPose=relPose;
    }
    
    //write the state of the reading and update all the particles using the motion model
    
    /*对于每一个粒子，都从里程计运动模型中采样，得到车子的初步估计位置  这一步对应于   里程计的更新 */

    for (ParticleVector::iterator it=m_particles.begin(); it!=m_particles.end(); it++){
      OrientedPoint& pose(it->pose);
      pose=m_motionModel.drawFromMotion(it->pose, relPose, m_odoPose);
    }

    // update the output file
    if (m_outputStream.is_open()){
      m_outputStream << setiosflags(ios::fixed) << setprecision(6);
      m_outputStream << "ODOM ";
      m_outputStream << setiosflags(ios::fixed) << setprecision(3) << m_odoPose.x << " " << m_odoPose.y << " ";
      m_outputStream << setiosflags(ios::fixed) << setprecision(6) << m_odoPose.theta << " ";
      m_outputStream << reading.getTime();
      m_outputStream << endl;
    }
    if (m_outputStream.is_open()){
      m_outputStream << setiosflags(ios::fixed) << setprecision(6);
      m_outputStream << "ODO_UPDATE "<< m_particles.size() << " ";
      for (ParticleVector::iterator it=m_particles.begin(); it!=m_particles.end(); it++){
	OrientedPoint& pose(it->pose);
	m_outputStream << setiosflags(ios::fixed) << setprecision(3) << pose.x << " " << pose.y << " ";
	m_outputStream << setiosflags(ios::fixed) << setprecision(6) << pose.theta << " " << it-> weight << " ";
      }
      m_outputStream << reading.getTime();
      m_outputStream << endl;
    }
    
    //invoke the callback
    /*回调函数  实际上什么都没做*/
    onOdometryUpdate();
    

    // accumulate the robot translation and rotation

	/*根据两次里程计的数据 计算出来机器人的线性位移和角度位移的累积值 
	   m_odoPose表示上一次的里程计位姿  relPose表示新的里程计的位姿*/
	
    OrientedPoint move=relPose-m_odoPose;
    move.theta=atan2(sin(move.theta), cos(move.theta));

	//统计机器人在进行激光雷达更新之前 走了多远的距离 以及　平移了多少的角度
    m_linearDistance+=sqrt(move*move);
    m_angularDistance+=fabs(move.theta);

	 /*
     * 如果机器人在走了m_distanceThresholdCheck这么远的距离都没有进行激光雷达的更新
     * 则需要进行报警。这个误差很可能是里程计或者激光雷达的BUG造成的。
     * 例如里程计数据出错 或者 激光雷达很久没有数据等等
     * 每次进行激光雷达的更新之后 m_linearDistance这个参数就会清零
     */
 
    // if the robot jumps throw a warning
    if (m_linearDistance>m_distanceThresholdCheck){
      cerr << "***********************************************************************" << endl;
      cerr << "********** Error: m_distanceThresholdCheck overridden!!!! *************" << endl;
      cerr << "m_distanceThresholdCheck=" << m_distanceThresholdCheck << endl;
      cerr << "Old Odometry Pose= " << m_odoPose.x << " " << m_odoPose.y 
	   << " " <<m_odoPose.theta << endl;
      cerr << "New Odometry Pose (reported from observation)= " << relPose.x << " " << relPose.y 
	   << " " <<relPose.theta << endl;
      cerr << "***********************************************************************" << endl;
      cerr << "** The Odometry has a big jump here. This is probably a bug in the   **" << endl;
      cerr << "** odometry/laser input. We continue now, but the result is probably **" << endl;
      cerr << "** crap or can lead to a core dump since the map doesn't fit.... C&G **" << endl;
      cerr << "***********************************************************************" << endl;
    }

	//更新 把当前的位置赋值给旧的位置
    m_odoPose=relPose;
    
    bool processed=false;

    // process a scan only if the robot has traveled a given distance or a certain amount of time has elapsed
  /*只有当机器人走过一定的距离  或者 旋转过一定的角度  或者过一段指定的时间才处理激光数据*/
	if (! m_count 
	|| m_linearDistance>=m_linearThresholdDistance 
	|| m_angularDistance>=m_angularThresholdDistance
    || (period_ >= 0.0 && (reading.getTime() - last_update_time_) > period_)){
      last_update_time_ = reading.getTime();      

      if (m_outputStream.is_open()){
	m_outputStream << setiosflags(ios::fixed) << setprecision(6);
	m_outputStream << "FRAME " <<  m_readingCount;
	m_outputStream << " " << m_linearDistance;
	m_outputStream << " " << m_angularDistance << endl;
      }
      
      if (m_infoStream)
	m_infoStream << "update frame " <<  m_readingCount << endl
		     << "update ld=" << m_linearDistance << " ad=" << m_angularDistance << endl;
      
      
      cerr << "Laser Pose= " << reading.getPose().x << " " << reading.getPose().y 
	   << " " << reading.getPose().theta << endl;
      
      
      //this is for converting the reading in a scan-matcher feedable form
       /*复制一帧数据 把激光数据转换为scan-match需要的格式*/
      assert(reading.size()==m_beams);
      double * plainReading = new double[m_beams];
      for(unsigned int i=0; i<m_beams; i++){
	plainReading[i]=reading[i];
      }
      m_infoStream << "m_count " << m_count << endl;
	  
	  //这个备份主要是用来储存的。

      RangeReading* reading_copy = 
              new RangeReading(reading.size(),
                               &(reading[0]),
                               static_cast<const RangeSensor*>(reading.getSensor()),
                               reading.getTime());

      if (m_count>0){

	  /*
            为每个粒子进行scanMatch，计算出来每个粒子的最优位姿，同时计算改最优位姿的得分和似然  对应于gmapping论文中的用最近的一次测量计算proposal的算法
            这里面除了进行scanMatch之外，还对粒子进行了权重的计算，并计算了粒子的有效区域 但不进行内存分配 内存分配在resample()函数中
            这个函数在gridslamprocessor.hxx里面。
       */ 
	     scanMatch(plainReading);
	   //至此 关于proposal的更新完毕了，接下来是计算权重
	   
	     if (m_outputStream.is_open()){
	         m_outputStream << "LASER_READING "<< reading.size() << " ";
	         m_outputStream << setiosflags(ios::fixed) << setprecision(2);
	         for (RangeReading::const_iterator b=reading.begin(); b!=reading.end(); b++){
	             m_outputStream << *b << " ";
	         }
	         OrientedPoint p=reading.getPose();
	          m_outputStream << setiosflags(ios::fixed) << setprecision(6);
	          m_outputStream << p.x << " " << p.y << " " << p.theta << " " << reading.getTime()<< endl;
	         m_outputStream << "SM_UPDATE "<< m_particles.size() << " ";
	         for (ParticleVector::const_iterator it=m_particles.begin(); it!=m_particles.end(); it++){
	           const OrientedPoint& pose=it->pose;
	            m_outputStream << setiosflags(ios::fixed) << setprecision(3) <<  pose.x << " " << pose.y << " ";
	            m_outputStream << setiosflags(ios::fixed) << setprecision(6) <<  pose.theta << " " << it-> weight << " ";
	        }
	     m_outputStream << endl;
	     }

		 //至此 关于proposal的更新完毕了，接下来是计算权重
	    onScanmatchUpdate();

		 /*
            由于scanMatch中对粒子的权重进行了更新，那么这个时候各个粒子的轨迹上的累计权重都需要重新计算
            这个函数即更新各个粒子的轨迹上的累计权重是更新
            GridSlamProcessor::updateTreeWeights(bool weightsAlreadyNormalized) 函数在gridslamprocessor_tree.cpp里面实现

        */
	    updateTreeWeights(false);
				
	    if (m_infoStream){
	       m_infoStream << "neff= " << m_neff  << endl;
	    }
	    if (m_outputStream.is_open()){
	       m_outputStream << setiosflags(ios::fixed) << setprecision(6);
	       m_outputStream << "NEFF " << m_neff << endl;
	    }

		/*
             * 粒子重采样  根据neff的大小来进行重采样  不但进行了重采样，也对地图进行更新
             * GridSlamProcessor::resample 函数在gridslamprocessor.hxx里面实现
             */

    	resample(plainReading, adaptParticles, reading_copy);
	
      } else {
	   m_infoStream << "Registering First Scan"<< endl;
	   for (ParticleVector::iterator it=m_particles.begin(); it!=m_particles.end(); it++){	
	       m_matcher.invalidateActiveArea();
	       m_matcher.computeActiveArea(it->map, it->pose, plainReading);
	       m_matcher.registerScan(it->map, it->pose, plainReading);
	  
	      // cyr: not needed anymore, particles refer to the root in the beginning!
	      
		   //为每个粒子创建路径的第一个节点。该节点的权重为0,父节点为it->node(这个时候为NULL)。
		  //因为第一个节点就是轨迹的根，所以没有父节点
		  
	      TNode* node=new	TNode(it->pose, 0., it->node,  0);
	      //node->reading=0;
          node->reading = reading_copy;
	      it->node=node;
	  
	  }
     }
      //		cerr  << "Tree: normalizing, resetting and propagating weights at the end..." ;
      updateTreeWeights(false);
      //		cerr  << ".done!" <<endl;
      
      delete [] plainReading;
      m_lastPartPose=m_odoPose; //update the past pose for the next iteration
      m_linearDistance=0;
      m_angularDistance=0;
      m_count++;
      processed=true;
      
      //keep ready for the next step
      for (ParticleVector::iterator it=m_particles.begin(); it!=m_particles.end(); it++){
	    it->previousPose=it->pose;
      }
      
    }
    if (m_outputStream.is_open())
      m_outputStream << flush;
    m_readingCount++;
    return processed;
  }
  
  
  std::ofstream& GridSlamProcessor::outputStream(){
    return m_outputStream;
  }
  
  std::ostream& GridSlamProcessor::infoStream(){
    return m_infoStream;
  }
  

  
  inline void GridSlamProcessor::scanMatch(const double* plainReading){
	/* sample a new pose from each scan in the reference */
   
	double sumScore=0;
	for (ParticleVector::iterator it=m_particles.begin(); it!=m_particles.end(); it++){
	  OrientedPoint corrected;
	  double score, l, s;
  /* 计算最优的粒子
  optimize 调用了 score 这个函数 （计算粒子得分）
  在score 函数里，首先计算障碍物的坐标phit，然后将phit转换成网格坐标iPhit
  计算光束上与障碍物相邻的非障碍物网格坐标pfree,pfrree由phit沿激光束方向移动一个网格的距离得到，将pfree转换成网格坐标ipfree（增量，并不是实际值）
  在iphit 及其附近8个（m_kernelSize:default=1）栅格（pr,对应自由栅格为pf）搜索最优可能是障碍物的栅格。
  最优准则： pr 大于某一阈值，pf小于该阈值，且pr栅格的phit的平均坐标与phit的距离bestMu最小。
  得分计算： s +=exp(-1.0/m_gaussianSigma*bestMu*besMu)  参考NDT算法,距离越大，分数越小，分数的较大值集中在距离最小值处，符合正态分布模型
  至此 score 函数结束并返回粒子（currentPose）得分，然后回到optimize函数
  optimize 干的事就是 currentPose 的位姿进行微调，前、后、左、右、左转、右转 共6次，然后选取得分最高的位姿，返回最终的得分
  */
	  score=m_matcher.optimize(corrected, it->map, it->pose, plainReading);
	  if (score>m_minimumScore){  //判断得分是否符合要求
		it->pose=corrected;
	  } else {
	  if (m_infoStream){
		m_infoStream << "Scan Matching Failed, using odometry. Likelihood=" << l <<std::endl;
		m_infoStream << "lp:" << m_lastPartPose.x << " "  << m_lastPartPose.y << " "<< m_lastPartPose.theta <<std::endl;
		m_infoStream << "op:" << m_odoPose.x << " " << m_odoPose.y << " "<< m_odoPose.theta <<std::endl;
	  }
	  }
  /*   likelihoodAndSocre 作用是计算粒子的权重和（l），如果出现匹配失败，则 l=noHit 	*/
	  m_matcher.likelihoodAndScore(s, l, it->map, it->pose, plainReading);
	  sumScore+=score;
	  it->weight+=l;
	  it->weightSum+=l;
   
  /* 计算可活动区域
	  //set up the selective copy of the active area
	  //by detaching the areas that will be updated
  computeActiveArea 用于计算每个粒子相应的位姿所扫描到的区域  
  计算过程首先考虑了每个粒子的扫描范围会不会超过子地图的大小，如果会，则resize地图的大小
  然后定义了一个activeArea 用于设置可活动区域，调用了gridLine() 函数,这个函数如何实现的，
  请参考百度文库那篇介绍。
  */
	  m_matcher.invalidateActiveArea();
	  m_matcher.computeActiveArea(it->map, it->pose, plainReading);
	}
	if (m_infoStream)
	  m_infoStream << "Average Scan Matching Score=" << sumScore/m_particles.size() << std::endl; 
  }
  

  
  int GridSlamProcessor::getBestParticleIndex() const{
    unsigned int bi=0;
    double bw=-std::numeric_limits<double>::max();
    for (unsigned int i=0; i<m_particles.size(); i++)
      if (bw<m_particles[i].weightSum){
	bw=m_particles[i].weightSum;
	bi=i;
      }
    return (int) bi;
  }

 /*
主要功能为归一化粒子的权重，同时计算出neff
*/
  inline void GridSlamProcessor::normalize()
  {
	//normalize the log m_weights
	double gain=1./(m_obsSigmaGain*m_particles.size());
	
	/*求所有粒子中的最大的权重*/
	double lmax= -std::numeric_limits<double>::max();
	for (ParticleVector::iterator it=m_particles.begin(); it!=m_particles.end(); it++)
	{
	  lmax=it->weight>lmax?it->weight:lmax;
	}
	//cout << "!!!!!!!!!!! maxwaight= "<< lmax << endl;
	
	/*权重以最大权重为中心的高斯分布*/
	m_weights.clear();
	double wcum=0;
	m_neff=0;
	for (std::vector<Particle>::iterator it=m_particles.begin(); it!=m_particles.end(); it++)
	{
	  m_weights.push_back(exp(gain*(it->weight-lmax)));
	  wcum+=m_weights.back();
	  //cout << "l=" << it->weight<< endl;
	}
	
	/*
	计算有效粒子数 和 归一化权重
	权重=wi/w
	neff = 1/w*w
	*/
	m_neff=0;
	for (std::vector<double>::iterator it=m_weights.begin(); it!=m_weights.end(); it++)
	{
	  *it=*it/wcum;
	  double w=*it;
	  m_neff+=w*w;
	}
	m_neff=1./m_neff;
	
  }

  
  double GridSlamProcessor::propagateWeights()
  {
	// don't calls this function directly, use updateTreeWeights(..) !
   
		  // all nodes must be resetted to zero and weights normalized
   
	  // the accumulated weight of the root
	  // 求所有根节点的累计权重之和
	  double lastNodeWeight=0;
   
	  // sum of the weights in the leafs
	  // 所有叶子节点的权重 也就是m_weights数组里面所有元素的和
	  double aw=0;
   
	  std::vector<double>::iterator w=m_weights.begin();
	  for (ParticleVector::iterator it=m_particles.begin(); it!=m_particles.end(); it++)
	  {
   
		  //求解所有叶子节点的累计权重
		  double weight=*w;
		  aw+=weight;
		  
		  //叶子节点的子节点累计权重就等于自己的权重 因为它没有子节点
		  //每一个粒子的路径都是从叶子节点开始的，得到了叶子节点，就得到了路径
		  TNode * n=it->node;
		  n->accWeight=weight;
   
		  lastNodeWeight+=propagateWeight(n->parent,n->accWeight);
		  
		  w++;
	  }
	  
	  if (fabs(aw-1.0) > 0.0001 || fabs(lastNodeWeight-1.0) > 0.0001)
	  {
		cerr << "ERROR: ";
		cerr << "root->accWeight=" << lastNodeWeight << "	 sum_leaf_weights=" << aw << endl;
		assert(0);		   
	  }
	  return lastNodeWeight;
  }
  
  /*
  @desc 粒子滤波器重采样。
  分为两步：
  1.需要重采样，则所有保留下来的粒子的轨迹都加上一个新的节点，然后进行地图更新。
  2.不需要冲采样，则所有的粒子的轨迹都加上一个新的节点，然后进行地图的更新
  在重采样完毕之后，会调用registerScan函数来更新地图
  */
  inline bool GridSlamProcessor::resample(const double* plainReading, int adaptSize, const RangeReading* reading)
  {
	
	bool hasResampled = false;
	
	/*备份老的粒子的轨迹  即保留叶子节点 在增加新节点的时候使用*/
	TNodeVector oldGeneration;
	for (unsigned int i=0; i<m_particles.size(); i++)
	{
	  oldGeneration.push_back(m_particles[i].node);
	}
	
	/*如果需要进行重采样*/
	if (m_neff<m_resampleThreshold*m_particles.size())
	{	  
	  
	  if (m_infoStream)
		m_infoStream  << "*************RESAMPLE***************" << std::endl;
	  
	  //采取重采样方法决定，哪些粒子会保留  保留的粒子会返回下标.里面的下标可能会重复，因为有些粒子会重复采样
	  //而另外的一些粒子会消失掉
	  uniform_resampler<double, double> resampler;
	  m_indexes=resampler.resampleIndexes(m_weights, adaptSize);
	  
	  if (m_outputStream.is_open())
	  {
		m_outputStream << "RESAMPLE "<< m_indexes.size() << " ";
		for (std::vector<unsigned int>::const_iterator it=m_indexes.begin(); it!=m_indexes.end(); it++)
		{
		  m_outputStream << *it <<	" ";
		}
		m_outputStream << std::endl;
	  }
	  
	  onResampleUpdate();
	  //BEGIN: BUILDING TREE
   
   
	  //重采样之后的粒子
	  ParticleVector temp;
	  unsigned int j=0;
	  
	  //要删除的粒子下标
	  std::vector<unsigned int> deletedParticles;		  //this is for deleteing the particles which have been resampled away.
	  
	  //枚举每一个要被保留的粒子
	  for (unsigned int i=0; i<m_indexes.size(); i++)
	  {
		//统计要被删除的粒子
		while(j<m_indexes[i])
		{
		  deletedParticles.push_back(j);
		  j++;
		}
		if (j==m_indexes[i])
		j++;
	
		//得到当前的保留的粒子
		Particle & p=m_particles[m_indexes[i]];
		
		//每一个需要保留下来的粒子都需要在路径中增加一个新的节点
		TNode* node=0;
		TNode* oldNode=oldGeneration[m_indexes[i]];
		
		//创建一个新的节点 改节点的父节点为oldNode
		node=new  TNode(p.pose, 0, oldNode, 0);
		node->reading=reading;
		
		//这个要保留下来的粒子，要保留的粒子的下标为m_indexs
		temp.push_back(p);
		temp.back().node=node;
		temp.back().previousIndex=m_indexes[i];
	  }
   
	  while(j<m_indexes.size())
	  {
		deletedParticles.push_back(j);
		j++;
	  }
	  
	  //把要删除的粒子的Node都删除掉，Node表示轨迹的起点(最新的点)
	  std::cerr <<	"Deleting Nodes:";
	  for (unsigned int i=0; i<deletedParticles.size(); i++)
	  {
		std::cerr <<" " << deletedParticles[i];
		delete m_particles[deletedParticles[i]].node;
		m_particles[deletedParticles[i]].node=0;
	  }
	  
	  std::cerr  << " Done" <<std::endl;
	  std::cerr << "Deleting old particles..." ;
	  std::cerr << "Done" << std::endl;
	  
	  //清楚全部的粒子 然后从tmp中读取保留下来的粒子
	  m_particles.clear();
	  
	  //枚举要保留下来的所有的粒子 每个粒子都需要更新地图
	  std::cerr << "Copying Particles and  Registering	scans...";
   
	  //对于保留下来的粒子进行更新 这里是可以并行化操作的。
	  //在并行化操作里面 m_particles.push_back()会报错 因此并行化 需要把push_back()提出来。
	  //在外面的的for循环进行
	  int tmp_size = temp.size();
  //#pragma omp parallel for
	  for(int i = 0; i<tmp_size;i++)
	  {
		  //对保留下来的粒子数据进行更新
		  //每个粒子的权重都设置为相同的值
		  temp[i].setWeight(0);
   
		  //为每个粒子更新running_scans
   
		  //增加了一帧激光数据 因此需要更新地图
		  m_matcher.registerScan(temp[i].map,temp[i].pose,plainReading);
		  //m_matcher.registerScan(temp[i].lowResolutionMap,temp[i].pose,plainReading);
	  }
   
	  //提取出来 防止并行优化时报错
	  for(int i = 0; i< tmp_size;i++)
		  m_particles.push_back(temp[i]);
   
	  std::cerr  << " Done" <<std::endl;
	  hasResampled = true;
	} 
	/*否则的话，进行扫描匹配*/
	else 
	{
	  //不进行重采样的话，权值不变。只为轨迹创建一个新的节点
   
	  //为每个粒子更新地图 同样可以并行化
	  int particle_size = m_particles.size();
  //#pragma omp parallel for
	  for(int i = 0; i < particle_size;i++)
	  {
		  //创建一个新的树节点
		  TNode* node = 0;
		  node = new TNode(m_particles[i].pose,0.0,oldGeneration[i],0);
   
		  //把这个节点接入到树中
		  node->reading = reading;
		  m_particles[i].node = node;
   
		  //更新各个例子的地图
		  m_matcher.invalidateActiveArea();
		  m_matcher.registerScan(m_particles[i].map, m_particles[i].pose, plainReading);
		  m_particles[i].previousIndex = i;
	  }
	  std::cerr<<std::endl;
   
   
	  std::cerr  << "Done" <<std::endl;
	  
	}
	//END: BUILDING TREE
	
	return hasResampled;
  }
 

  void GridSlamProcessor::onScanmatchUpdate(){}
  void GridSlamProcessor::onResampleUpdate(){}
  void GridSlamProcessor::onOdometryUpdate(){}

  
};// end namespace




