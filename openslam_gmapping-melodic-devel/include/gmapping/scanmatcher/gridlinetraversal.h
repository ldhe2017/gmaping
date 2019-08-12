#ifndef GRIDLINETRAVERSAL_H
#define GRIDLINETRAVERSAL_H

#include <cstdlib>
#include <gmapping/utils/point.h>

namespace GMapping {

typedef struct {
  int     num_points;
  IntPoint*  points;
} GridLineTraversalLine;

struct GridLineTraversal {
  inline static void gridLine( IntPoint start, IntPoint end, GridLineTraversalLine *line ) ;
  inline static void gridLineCore( IntPoint start, IntPoint end, GridLineTraversalLine *line ) ;

};
//Bresenham算法
void GridLineTraversal::gridLineCore( IntPoint start, IntPoint end, GridLineTraversalLine *line )
{
  int   dx, dy;            //横纵坐标间距
  int	incr1, incr2;      //p_m增量
  int   d ;                 //p_m
  int x, y, xend, yend;    //直线增长的首末端点坐标
  int xdirflag, ydirflag;   //横纵坐标增长方向
  int cnt = 0;              //直线过点的点的序号

  dx = abs(end.x-start.x); dy = abs(end.y-start.y);
  
  if (dy <= dx) {                     //斜率K的绝对值|K|<1时，在x方向进行单位步进
    d = 2*dy - dx;                   //初始点p_m0值
	incr1 = 2 * dy;                  //情况1
    incr2 = 2 * (dy - dx);           //情况2
    if (start.x > end.x) {           //起点横坐标比终点横坐标大，xdirflag = -1（负号可以理解为增长方向与直线始终点方向相反）
      x = end.x; y = end.y;          //设置增长起点，注意这里要选择小横坐标作为起点，用于增量时好理解
      ydirflag = (-1);               //此时默认(希望)纵坐标也是start.y > end.y的情况
      xend = start.x;
    } else {
      x = start.x; y = start.y;
      ydirflag = 1;
      xend = end.x;
    }
    line->points[cnt].x=x;               //加入起点坐标
    line->points[cnt].y=y;
    cnt++;
    if (((end.y - start.y) * ydirflag) > 0) {    //是预料到的情况，start.y > end.y
      while (x < xend) {                          //开始进行增量递增

	  //x > xend，y > yend，处理向右上角爬升的直线群，下图1号所示直线
	  
	x++;
	if (d <0) {
	  d+=incr1;
	} else {
	  y++; d+=incr2;                 //纵坐标向正方向增长
	}
	line->points[cnt].x=x;          //添加新的点
	line->points[cnt].y=y;
	cnt++;
      }
    } else {                        //x > xend，y < yend，处理向右下角降落的直线群，下图2号所示直线
      while (x < xend) {
	x++;
	if (d <0) {
	  d+=incr1;
	} else {
	  y--; d+=incr2;              //纵坐标向负方向增长
	}
	line->points[cnt].x=x;
	line->points[cnt].y=y;
	cnt++;
      }
    }		
  } else {             //dy > dx，当斜率k的绝对值|k|>1时，在y方向进行单位步进
    d = 2*dx - dy;                        //P_m0初值
    incr1 = 2*dx; incr2 = 2 * (dx - dy);  //形同算法推导情况1
    if (start.y > end.y) {
      y = end.y; x = end.x;               //取最小的纵坐标作为起点
      yend = start.y;
      xdirflag = (-1);                   //期望start.x > end.x
    } else {
      y = start.y; x = start.x;
      yend = end.y;
      xdirflag = 1;
    }
    line->points[cnt].x=x;               //添加起点
    line->points[cnt].y=y;
    cnt++;
    if (((end.x - start.x) * xdirflag) > 0) {     //x > xend ，y > yend，处理向右上角爬升的直线群，下图3号所示直线
      while (y < yend) {
	y++;
	if (d <0) {
	  d+=incr1;
	} else {
	  x++; d+=incr2;                 //横坐标向正方向增长
	}
	line->points[cnt].x=x;            //添加新的点
	line->points[cnt].y=y;
	cnt++;
      }
    } else {                          //x < xend，y > yend，处理向左上角爬升的直线群，下图4号所示直线
      while (y < yend) {
	y++;
	if (d <0) {
	  d+=incr1;
	} else { 
	  x--; d+=incr2;                 //横坐标向负方向增长
	}
	line->points[cnt].x=x;
	line->points[cnt].y=y;
	cnt++;
      }
    }
  }
  line->num_points = cnt;             //记录添加所有点的数目
}

void GridLineTraversal::gridLine( IntPoint start, IntPoint end, GridLineTraversalLine *line ) {
  int i,j;
  int half;
  IntPoint v;
  gridLineCore( start, end, line );
  if ( start.x!=line->points[0].x ||
       start.y!=line->points[0].y ) {
    half = line->num_points/2;
    for (i=0,j=line->num_points - 1;i<half; i++,j--) {
      v = line->points[i];
      line->points[i] = line->points[j];
      line->points[j] = v;
    }
  }
}

};

#endif
