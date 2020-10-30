/* The copyright in this software is being made available under the BSD
 * Licence, included below.  This software may be subject to other third
 * party and contributor rights, including patent rights, and no such
 * rights are granted under this licence.
 *
 * Copyright (c) 2017-2018, ISO/IEC
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the distribution.
 *
 * * Neither the name of the ISO/IEC nor the names of its contributors
 *   may be used to endorse or promote products derived from this
 *   software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include "RAHT.h"

#include <cassert>
#include <cinttypes>
#include <climits>
#include <cstddef>
#include <utility>
#include <vector>
#include <stdio.h>

#include "PCCTMC3Common.h"
#include "PCCMisc.h"

namespace pcc {

//============================================================================

struct UrahtNode {
  int64_t pos;
  int weight;
  Qps qp;
};

//============================================================================
// remove any non-unique leaves from a level in the uraht tree

int
reduceUnique(
  int numNodes,
  int numAttrs,
  std::vector<UrahtNode>* weightsIn,
  std::vector<UrahtNode>* weightsOut,
  std::vector<int>* attrsIn,
  std::vector<int>* attrsOut)
{
  // process a single level of the tree
 //   numnodes:当前层中所包含的点数
 //   numAttrs:每个点包含的属性信息的个数，例如rgb即为3个信息
 //   weightsIn:点云序列中所有的点的包含的位置信息、权重信息、量化参数信息
 //   weightsOut:构建当前层的点的信息
 //   attrsIn:点云序列中所有的点的属性信息
	//attrsOut:当前层的属性信息

  int64_t posPrev = -1; //设置一个初始坐标为-1
  auto weightsInWrIt = weightsIn->begin();//指向点云序列中的点的首地址
  auto weightsInRdIt = weightsIn->cbegin();//指向点云序列中的点的首地址，cbegin不能修改指针指向的元素
  auto attrsInWrIt = attrsIn->begin();//指向点云属性信息的首地址
  auto attrsInRdIt = attrsIn->begin();//指向点云属性信息的首地址
  for (int i = 0; i < numNodes; i++) {
    const auto& node = *weightsInRdIt++;//通过node记录当前节点的信息

    // copy across unique nodes
    if (node.pos != posPrev) {//当前点没有重复出现：即莫顿码信息与上一个点不同
      posPrev = node.pos;//记录当前点的莫顿码信息
      *weightsInWrIt++ = node;//输出当前点的信息
      for (int k = 0; k < numAttrs; k++)
        *attrsInWrIt++ = *attrsInRdIt++;//输出当前点的属性信息
      continue;
    }

    // duplicate node
	//当前点重复出现，即当前点与上一个点的莫顿码信息相同，此时需进行去除重复点的操作
    (weightsInWrIt - 1)->weight += node.weight;//更新上一个点的权重信息=当前点的权重信息+上一个点的权重信息
    weightsOut->push_back(node);//在weightsOut中输出当前点的信息
    for (int k = 0; k < numAttrs; k++) {
      *(attrsInWrIt - numAttrs + k) += *attrsInRdIt;//更新上一个点的属性信息值=当前点的属性值+上一个点对应的属性值
      attrsOut->push_back(*attrsInRdIt++);//在attrsOut中输出当前点的属性信息值
    }
  }

  // number of nodes in next level
  return std::distance(weightsIn->begin(), weightsInWrIt);
}

//============================================================================
// Split a level of values into sum and difference pairs.

int
reduceLevel(
  int level,
  int numNodes,
  int numAttrs,
  std::vector<UrahtNode>* weightsIn,
  std::vector<UrahtNode>* weightsOut,
  std::vector<int>* attrsIn,
  std::vector<int>* attrsOut)
{
  // process a single level of the tree
	// 通过给当前点和上一点的莫顿码右移相同的位数，判断这两个节点是否属于同一父节点，从而进行点的合并
  //   level:当前层的层数
  //   numnodes:上一层中所包含的点数
  //   numAttrs:每个点包含的属性信息的个数，例如rgb即为3个信息
  //   weightsIn:上一层中的位置信息，经过构建树之后更新为当前层中所有的点的包含的位置信息、权重信息、量化参数信息
  //   weightsOut:同一父节点同行/列的另一个点的信息
  //   attrsIn:点云序列上一层点的属性信息，更新为当前层的点的属性信息
  //   attrsOut:同一父节点同行/列的另一个点的属性信息
  int64_t posPrev = -1;//设置一个初始坐标为-1
  auto weightsInWrIt = weightsIn->begin();//指向点云序列中的点的首地址
  auto weightsInRdIt = weightsIn->cbegin();//指向点云序列中的点的首地址，cbegin不能修改指针指向的元素
  auto attrsInWrIt = attrsIn->begin();//指向点云属性信息的首地址
  auto attrsInRdIt = attrsIn->begin();//指向点云属性信息的首地址
  for (int i = 0; i < numNodes; i++) {//按照莫顿码顺序遍历点云中的点
    auto& node = *weightsInRdIt++;//node记录当前点的信息
    bool newPair = (posPrev ^ node.pos) >> level != 0;
	//判断当前点和莫顿码比它小的前一个点是否是邻居，newPair为真时不是邻居，为假时是邻居
    posPrev = node.pos;
    posPrev = node.pos;//posPrev记录当前点的莫顿码信息
    //当前点与前一个点不为邻居关系时，将当前点几何信息存入weightsLf，属性信息存入attrsLf
    if (newPair) {
      *weightsInWrIt++ = node;//输出当前点的信息在weightsIn中
      for (int k = 0; k < numAttrs; k++)
        *attrsInWrIt++ = *attrsInRdIt++;//输出当前点的属性信息在attrsIn中
    } else {
      //是邻居时，将这两个节点的权重相加存入weightsLf，属性相加存入attrsLf，莫顿码为前一个点的莫顿码值。将当前点几何信息存入weightsHf,属性信息存入attrsHf
      auto& left = *(weightsInWrIt - 1);//记录当前点的上一个节点的信息
      left.weight += node.weight;//更新上一个节点的权重信息=上一个节点的权重信息+当前节点的权重信息
      left.qp[0] = (left.qp[0] + node.qp[0]) >> 1;//更新上一个节点的量化参数信息=当前点+上一个节点的权重信息/2
      left.qp[1] = (left.qp[1] + node.qp[1]) >> 1;//更新上一个节点的量化参数信息=当前点+上一个节点的权重信息/2
      weightsOut->push_back(node);//输出当前点到weightOut中

      for (int k = 0; k < numAttrs; k++) {
        *(attrsInWrIt - numAttrs + k) += *attrsInRdIt;//更新上一个点的属性信息值=当前点的属性值+上一个点对应的属性值
        attrsOut->push_back(*attrsInRdIt++);//输出当前点的属性值到attrsOut中
      }
    }
  }

  // number of nodes in next level 返回该层没有邻居的点数
  return std::distance(weightsIn->begin(), weightsInWrIt);
}

//============================================================================
// Merge sum and difference values to form a tree.

void
expandLevel(
  int level,
  int numNodes,
  int numAttrs,
  std::vector<UrahtNode>* weightsIn,   // expand by numNodes before expand
  std::vector<UrahtNode>* weightsOut,  // shrink after expand
  std::vector<int>* attrsIn,
  std::vector<int>* attrsOut)
{
  // 通过给当前点和上一点的莫顿码右移相同的位数，判断这两个节点是否属于同一父节点，从而进行点的合并
  //   level:当前层的层数
  //   numnodes:当前层中所包含的重复点个数 当前层的点数=上一层中点的个数+当前层中重复点的个数
  //   numAttrs:每个点包含的属性信息的个数，例如rgb即为3个信息
  //   weightsIn:当前层中所有的点的包含的位置信息、权重信息、量化参数信息
  //   weightsOut:当前层中重复点的信息
  //   attrsIn:当前层的点的属性信息
  //   attrsOut:当前层中重复点的属性信息
  if (numNodes == 0)//判断该层中点数是否为0，若为0，则直接退出该函数
    return;

  // process a single level of the tree
  auto weightsInWrIt = weightsIn->rbegin();//读取当前层最后一个点的信息
  auto weightsInRdIt = std::next(weightsIn->crbegin(), numNodes);//读取上层最后一个点的信息
  auto weightsOutRdIt = weightsOut->crbegin();//所有重复点中最后一个重复点的信息
  auto attrsInWrIt = attrsIn->rbegin();//当前层最后一个点的最后一个属性信息
  auto attrsInRdIt = std::next(attrsIn->crbegin(), numNodes * numAttrs);//层中最后一个点的第一个属性信息值
  auto attrsOutRdIt = attrsOut->crbegin();//重复点中最后一个重复点的最后一个属性信息值
  for (int i = 0; i < numNodes;) {
    //判断当前点与weightsHf中从后开始遍历的点是否是邻居关系，isPair为真表示是邻居关系，为假表示不是邻居关系
    bool isPair = (weightsOutRdIt->pos ^ weightsInRdIt->pos) >> level == 0;
    if (!isPair) {//如果不属于邻居节点
      *weightsInWrIt++ = *weightsInRdIt++;//将当前层中上一层点的信息赋给当前重复点
      for (int k = 0; k < numAttrs; k++)
        *attrsInWrIt++ = *attrsInRdIt++;//将当前层中上一层点的属性信息赋当前的重复点
      continue;
    }

	//属于邻居节点，则当前点的属性减去attrsHf从后开始遍历点的属性值得到当前点邻居节点的属性值
    // going to process a pair
    i++;//重复点计数+1

    // Out node is inserted before In node.
    const auto& nodeDelta = *weightsInWrIt++ = *weightsOutRdIt++;//将当前重复点的信息赋给当前层的点
    auto curAttrIt = attrsInWrIt;//读取当前层中点的属性值
    for (int k = 0; k < numAttrs; k++)
      *attrsInWrIt++ = *attrsOutRdIt++;//重复点的值赋给当前点

    // move In node to correct position, subtracting delta
    *weightsInWrIt = *weightsInRdIt++;//上一层中对应点的位置
    (weightsInWrIt++)->weight -= nodeDelta.weight;//求出上一层中点在当前层的权重
    for (int k = 0; k < numAttrs; k++) {
      *attrsInWrIt = *attrsInRdIt++;//上一层中对应点的属性信息
      *attrsInWrIt++ -= *curAttrIt++;//求出上一层中点在当前层的属性信息
    }
  }
}

//============================================================================
// Search for neighbour with @value in the ordered list [first, last).
//
// If distance is positive, search [from, from+distance].
// If distance is negative, search [from-distance, from].

template<typename It, typename T, typename T2, typename Cmp>
It
findNeighbour(It first, It last, It from, T value, T2 distance, Cmp compare)
{
  It start = first;
  It end = last;

  if (distance >= 0) {
    start = from;
    if ((distance + 1) < std::distance(from, last))
      end = std::next(from, distance + 1);
  } else {
    end = from;
    if ((-distance) < std::distance(first, from))
      start = std::prev(from, -distance);
  }

  auto found = std::lower_bound(start, end, value, compare);
  if (found == end)
    return last;
  return found;
}

//============================================================================
// Find the neighbours of the node indicated by @t between @first and @last.
// The position weight of each found neighbour is stored in two arrays.

template<typename It>
void
findNeighbours(
  It first,
  It last,
  It it,
  int level,
  uint8_t occupancy,
  int parentNeighIdx[19],
  int parentNeighWeights[19])
{
	//1-19分别标记固定位置的邻居节点，[1]表示它本身，共面，共线的邻居
  //对于当前节点，其邻居有6个共面，12个共线包括当前父节点本身，一共最多有19个邻居，neighMasks表示当前父节点中
  //与邻居共面或共线的有效子节点，1表示有效，0表示无效。比如15（00001111），表示对于左面邻居，当前父节点内索引为0,1,2,3的子节点是有效子节点。
  static const uint8_t neighMasks[19] = {255, 15, 240, 51, 204, 85,  170,
                                         3,   12, 5,   10, 48,  192, 80,
                                         160, 17, 34,  68, 136};
  //用8bit表示当前节点可预测的子节点的位置，8bit分别表示当前父块中子节点的位置信息，1表示可进行预测，0表示不能预测该子块

  // current position (discard extra precision)
  int64_t cur_pos = it->pos >> level;//计算当前父节点的地址

  // the position of the parent, offset by (-1,-1,-1)
  int64_t base_pos = morton3dAdd(cur_pos, -1ll);//对当前父节点坐标进行偏移

  // these neighbour offsets are relative to base_pos
  static const uint8_t neighOffset[19] = {0,  3,  35, 5,  21, 6, 14, 1,  17, 2,
                                          10, 33, 49, 34, 42, 4, 12, 20, 28};

   //表示每一个邻居节点相对于当前节点的几何坐标偏移量
  // special case for the direct parent (no need to search);
  parentNeighIdx[0] = std::distance(first, it);
  parentNeighWeights[0] = it->weight;//第0个邻居为当前父节点本身

  for (int i = 1; i < 19; i++) {
    // Only look for neighbours that have an effect
    if (!(occupancy & neighMasks[i])) {
      parentNeighIdx[i] = -1;
      continue;
    }
    // 仅当父节点内对邻居有效的子节点占据时才会搜索该邻居，比如当索引为0, 1, 2,3的子节点为空时，不会搜索左侧共面邻居

    // compute neighbour address to look for
    // the delta between it and the current position is
    int64_t neigh_pos = morton3dAdd(base_pos, neighOffset[i]);//根据邻居节点对父节点的坐标偏移量计算邻居节点的莫顿码值
    int64_t delta = neigh_pos - cur_pos;//计算邻居节点相对于当前父节点的位置残差


    // find neighbour
    auto found = findNeighbour(
      first, last, it, neigh_pos, delta,
      [=](decltype(*it)& candidate, int64_t neigh_pos) {
        return (candidate.pos >> level) < neigh_pos;
      });//根据邻居与当前父节点的距离残差搜索邻居

    if (found == last) {
      parentNeighIdx[i] = -1;
      continue;
    }//当找到的邻居地址不存在时，该邻居不存在，索引赋为-1

    if ((found->pos >> level) != neigh_pos) {
      parentNeighIdx[i] = -1;
      continue;
    }//当根据距离残差找到的邻居莫顿码右移level以后不等于邻居地址，该邻居不存在，索引赋为-1

    parentNeighIdx[i] = std::distance(first, found);//计算邻居索引
    parentNeighWeights[i] = found->weight;//计算邻居权重
  }//计算其余18个邻居
}

//============================================================================
// Generate the spatial prediction of a block.

template<typename It>
void
intraDcPred(
  int numAttrs,
  const int neighIdx[19],
  const int neighWeights[19],
  int occupancy,
  It first,
  FixedPoint predBuf[][8])
{
  //每一个邻居只对当前父节点内与其共面或共线的子节点进行属性预测，当predMasks的值转为二进制以后，为1的子节点可用该邻居进行预测，为0的子节点不进行预测，
  //比如，predMasks[1]=15，即为00001111，表示左侧共面邻居只用来预测索引为0,1,2,3的子节点，而predMasks[0]=255表示当前父节点预测所有子节点；
  static const uint8_t predMasks[19] = {255, 15, 240, 51, 204, 85,  170,
                                        3,   12, 5,   10, 48,  192, 80,
                                        160, 17, 34,  68, 136};
  //预测权重
  //预测权重，当前父节点对子节点的预测权重为4，共面邻居预测权重为2，共线邻居预测权重为1
  static const int predWeight[19] = {4, 2, 2, 2, 2, 2, 2, 1, 1, 1,
                                     1, 1, 1, 1, 1, 1, 1, 1, 1};

    //替换除法操作的查询表
  static const int kDivisors[25] = {8192, 6554, 5461, 4681, 4096, 3641, 3277,
                                    2979, 2731, 2521, 2341, 2185, 2048, 1928,
                                    1820, 1725, 1638, 1560, 1489, 1425, 1365,
                                    1311, 1260, 1214, 1170};

  int weightSum[8] = {-4, -4, -4, -4, -4, -4, -4, -4};//预测权重之和

  std::fill_n(&predBuf[0][0], 8 * numAttrs, FixedPoint(0));

  int64_t neighValue[3];
  //去除无效邻居的阈值
  int64_t limitLow = 0;//预测阈值下限
  int64_t limitHigh = 0;//预测阈值上限
  for (int i = 0; i < 19; i++) {
    if (neighIdx[i] == -1)//该点不存在父邻居节点
      continue;//退出该次循环

    auto neighValueIt = std::next(first, numAttrs * neighIdx[i]);
    for (int k = 0; k < numAttrs; k++)
      neighValue[k] = *neighValueIt++;

    // skip neighbours that are outside of threshold limits
	//判断邻居父节点是否在阈值内，若不在，则该点不能作为预测值去预测子块的信息
    //当邻居与当前父节点的属性值比例小于limitLow或大于limitHigh时，去除该邻居，这里只考虑Y分量的属性值
    if (i) {
      if (10 * neighValue[0] <= limitLow || 10 * neighValue[0] >= limitHigh)
        continue;
    } else {
      constexpr int ratioThreshold1 = 2;
      constexpr int ratioThreshold2 = 25;
      limitLow = ratioThreshold1 * neighValue[0];//计算邻居属性的最小阈值
      limitHigh = ratioThreshold2 * neighValue[0];//计算邻居属性的最大阈值
    }

    // apply weighted neighbour value to masked positions
    for (int k = 0; k < numAttrs; k++)
      neighValue[k] *= predWeight[i] << pcc::FixedPoint::kFracBits;

    //加权求预测属性值
    int mask = predMasks[i] & occupancy;
    for (int j = 0; mask; j++, mask >>= 1) {
      if (mask & 1) {
        weightSum[j] += predWeight[i];
        for (int k = 0; k < numAttrs; k++)
          predBuf[k][j].val += neighValue[k];
      }
    }
  }

  // normalise
  //对预测属性值进行归一化
  FixedPoint div;
  for (int i = 0; i < 8; i++, occupancy >>= 1) {
    if (occupancy & 1) {
      div.val = kDivisors[weightSum[i]];
      for (int k = 0; k < numAttrs; k++)
        predBuf[k][i] *= div;
    }
  }
}

//============================================================================
// Encapsulation of a RAHT transform stage.

class RahtKernel {
public:
  RahtKernel(int weightLeft, int weightRight)
  {
    uint64_t w = weightLeft + weightRight;
    uint64_t isqrtW = irsqrt(w);
    _a.val =
      (isqrt(uint64_t(weightLeft) << (2 * _a.kFracBits)) * isqrtW) >> 40;
    _b.val =
      (isqrt(uint64_t(weightRight) << (2 * _b.kFracBits)) * isqrtW) >> 40;
  }

  void fwdTransform(
    FixedPoint left, FixedPoint right, FixedPoint* lf, FixedPoint* hf)
  {
    FixedPoint a = _a, b = _b;
    // lf = left * a + right * b
    // hf = right * a - left * b

    *lf = right;
    *lf *= b;
    *hf = right;
    *hf *= a;

    a *= left;
    b *= left;

    *lf += a;
    *hf -= b;
  }

  void invTransform(
    FixedPoint lf, FixedPoint hf, FixedPoint* left, FixedPoint* right)
  {
    FixedPoint a = _a, b = _b;

    *left = lf;
    *left *= a;
    *right = lf;
    *right *= b;

    b *= hf;
    a *= hf;

    *left -= b;
    *right += a;
  }

private:
  FixedPoint _a, _b;
};

//============================================================================
// In-place transform a set of sparse 2x2x2 blocks each using the same weights
//RAHT变换
template<class Kernel>
void
fwdTransformBlock222(
  int numBufs, FixedPoint buf[][8], int weights[8 + 8 + 8 + 8])
{
  static const int a[4 + 4 + 4] = {0, 2, 4, 6, 0, 4, 1, 5, 0, 1, 2, 3};//RAHT变换顺序
  static const int b[4 + 4 + 4] = {1, 3, 5, 7, 2, 6, 3, 7, 4, 5, 6, 7};
  for (int i = 0, iw = 0; i < 12; i++, iw += 2) {
    int i0 = a[i];//输出进行变换的块的索引
    int i1 = b[i];

    if (weights[iw] + weights[iw + 1] == 0)//若两个节点的权重之和为0，则该点不进行RAHT变换
      continue;

    // only one occupied, propagate to next level
	//如果进行RAHT变换的两个节点，只有一个节点被占据，那么直接交换两个节点的位置，传入下一层，准备进行下次变换
    if (!weights[iw] || !weights[iw + 1]) {
      if (!weights[iw]) {
        for (int k = 0; k < numBufs; k++)
          std::swap(buf[k][i0], buf[k][i1]);
      }
      continue;
    }

	//两个节点都被占据，进行RAHT变换
    // actual transform
    Kernel kernel(weights[iw], weights[iw + 1]);
    for (int k = 0; k < numBufs; k++) {
      auto& bufk = buf[k];
      kernel.fwdTransform(bufk[i0], bufk[i1], &bufk[i0], &bufk[i1]);
    }
  }
}

//============================================================================
// In-place inverse transform a set of sparse 2x2x2 blocks each using the
// same weights

template<class Kernel>
void
invTransformBlock222(
  int numBufs, FixedPoint buf[][8], int weights[8 + 8 + 8 + 8])
{
  static const int a[4 + 4 + 4] = {0, 2, 4, 6, 0, 4, 1, 5, 0, 1, 2, 3};
  static const int b[4 + 4 + 4] = {1, 3, 5, 7, 2, 6, 3, 7, 4, 5, 6, 7};
  for (int i = 11, iw = 22; i >= 0; i--, iw -= 2) {
    int i0 = a[i];
    int i1 = b[i];

    if (weights[iw] + weights[iw + 1] == 0)
      continue;

    // only one occupied, propagate to next level
    if (!weights[iw] || !weights[iw + 1]) {
      if (!weights[iw]) {
        for (int k = 0; k < numBufs; k++)
          std::swap(buf[k][i0], buf[k][i1]);
      }
      continue;
    }

    // actual transform
    Kernel kernel(weights[iw], weights[iw + 1]);
    for (int k = 0; k < numBufs; k++) {
      auto& bufk = buf[k];
      kernel.invTransform(bufk[i0], bufk[i1], &bufk[i0], &bufk[i1]);
    }
  }
}

//============================================================================
// expand a set of eight weights into three levels

void
mkWeightTree(int weights[8 + 8 + 8 + 8])
{
  int* in = &weights[0];
  int* out = &weights[8];

  for (int i = 0; i < 4; i++) {
    out[0] = out[4] = in[0] + in[1];
    if (!in[0] || !in[1])
      out[4] = 0;  // single node, no high frequencies
    in += 2;
    out++;
  }
  out += 4;
  for (int i = 0; i < 4; i++) {
    out[0] = out[4] = in[0] + in[1];
    if (!in[0] || !in[1])
      out[4] = 0;  // single node, no high frequencies
    in += 2;
    out++;
  }
  out += 4;
  for (int i = 0; i < 4; i++) {
    out[0] = out[4] = in[0] + in[1];
    if (!in[0] || !in[1])
      out[4] = 0;  // single node, no high frequencies
    in += 2;
    out++;
  }
}

//============================================================================
// Invoke mapFn(coefIdx) for each present coefficient in the transform

template<class T>
void
scanBlock(int weights[8 + 8 + 8 + 8], T mapFn)
{
  static const int8_t kRahtScanOrder[] = {0, 4, 2, 1, 6, 5, 3, 7};

  // there is always the DC coefficient (empty blocks are not transformed)
  mapFn(0);

  for (int i = 1; i < 8; i++) {
    if (!weights[24 + kRahtScanOrder[i]])
      continue;

    mapFn(kRahtScanOrder[i]);
  }
}

//============================================================================
// Tests if two positions are siblings at the given tree level

static bool
isSibling(int64_t pos0, int64_t pos1, int level)
{
  return ((pos0 ^ pos1) >> level) == 0;
}

//============================================================================
// Core transform process (for encoder/decoder)

template<bool isEncoder>
void
uraht_process(
  bool raht_prediction_enabled_flag,
  const int predictionThreshold[2],
  const QpSet& qpset,
  const Qps* pointQpOffsets,
  int numPoints,
  int numAttrs,
  int64_t* positions,
  int* attributes,
  int32_t* coeffBufIt)
{
  //raht_prediction_enabled_flag:属性预测开关
  //predictionThreshold[2]:属性预测阈值，判断父节点个数和祖父节点个数
  //qpset:量化参数
  //pointQpOffsets:量化偏移参数
  //numPoints:点云序列中的点数
  //numAttrs:属性信息个数
  //positions:位置，即莫顿码信息
  //attributes:属性信息
  //coeffBufIt:变换系数信息
  // coefficients are stored in three planar arrays.  coeffBufItK is a set
  // of iterators to each array.
  //将三个属性信息对应的变换系数的数组的首地址写入coeffBufItK[3]中
  int32_t* coeffBufItK[3] = {
    coeffBufIt,
    coeffBufIt + numPoints,
    coeffBufIt + numPoints * 2,
  };

  
  if (numPoints == 1) {//若当前点云序列中的点数为1
    auto quantizers = qpset.quantizers(0, pointQpOffsets[0]);//选择量化参数，Y分量和Cb、Cr之间可能存在量化偏移
    for (int k = 0; k < numAttrs; k++) {
      auto& q = quantizers[std::min(k, int(quantizers.size()) - 1)];//选择量化参数，Y分量选择quantizers[0],Cb、Cr分量选择quantizers[1];

      if (isEncoder) {//编码
        auto coeff = attributes[k];//把当前点的属性信息值赋给coeff
        assert(coeff <= INT_MAX && coeff >= INT_MIN);//判断当前属性信息值在[-2^32,2^32]之间
        *coeffBufItK[k]++ = coeff =
          q.quantize(coeff << kFixedPointAttributeShift);//对当前属性信息值coeff进行量化，并存入coeffBufItk内
        attributes[k] =
          divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift);//对量化之后的属性值进行反量化并写入attributes中
      } else {//解码
        int64_t coeff = *coeffBufItK[k]++;//读取量化之后的信息
        attributes[k] =
          divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift);//进行反量化
      }
    }
    return;
  }

  std::vector<UrahtNode> weightsLf, weightsHf;
  //创建RAHT变换树的左右树枝，数据类型为UrahtNode类型，包含了位置信息：莫顿码，权重信息，量化参数信息
  //在构建RAHT变换树时，左侧的树不断更新状态信息，右侧的树记录重复点的信息
  std::vector<int> attrsLf, attrsHf;//创建存放属性信息的RAHT变换树左右树枝

  weightsLf.reserve(numPoints);//设置树枝的大小=点云序列的点数
  attrsLf.reserve(numPoints * numAttrs);//树枝属性树枝的大小=点云序列的点数*每个点包含的属性值的个数

  int regionQpShift = 4;//设置原始的量化偏移量为4

  // copy positions into internal form
  // todo(df): lift to api
  for (int i = 0; i < numPoints; i++) {//在左侧树枝中按莫顿序存放莫顿码信息，权重信息，量化偏移量信息，每个点的初始权重设置为1
    weightsLf.emplace_back(UrahtNode{positions[i],
                                     1,
                                     {pointQpOffsets[i][0] << regionQpShift,
                                      pointQpOffsets[i][1] << regionQpShift}});
    for (int k = 0; k < numAttrs; k++) {
      attrsLf.push_back(attributes[i * numAttrs + k]);//在RAHT属性信息树左侧树枝依次存放每个点的各个属性信息
    }
  }

  weightsHf.reserve(numPoints);//设置树枝的大小=点云序列的点数
  attrsHf.reserve(numPoints * numAttrs);//树枝属性树枝的大小=点云序列的点数*每个点包含的属性值的个数

  // ascend tree
  std::vector<int> levelHfPos;//创建levelHfPos存放每层树中包含的点的个数
  //构建RAHT预测树

  for (int level = 0, numNodes = weightsLf.size(); numNodes > 1; level++) {
    levelHfPos.push_back(weightsHf.size());
    //存储每一层与前一个点是邻居的成对的点数
    if (level == 0) {//若当前层为第0层
      // process any duplicate points在第一层去除原始点云中的重复点
	//使用函数reduceUnique构建第0层的树，返回该层中包含的点的个数
      numNodes = reduceUnique(
        numNodes, numAttrs, &weightsLf, &weightsHf, &attrsLf, &attrsHf);
    } else {
      // normal level reduction 判断每一层点的邻居关系
		//逐层构建用于RAHT变换的树
      numNodes = reduceLevel(
        level, numNodes, numAttrs, &weightsLf, &weightsHf, &attrsLf, &attrsHf);
    }
  }

  assert(weightsLf[0].weight == numPoints);//保证点云中所有的点最后都会合并在第一个点中，即它的权重等于所有点的权重之和

  // reconstruction buffers
  std::vector<int> attrRec, attrRecParent;//创建容器存放重建信息
  attrRec.resize(numPoints * numAttrs);//设置容器大小
  attrRecParent.resize(numPoints * numAttrs);

  std::vector<int> attrRecUs, attrRecParentUs;//创建容器存放重建信息
  attrRecUs.resize(numPoints * numAttrs);//设置容器大小
  attrRecParentUs.resize(numPoints * numAttrs);

  std::vector<UrahtNode> weightsParent;//创建容器存放父节点的权重信息
  weightsParent.reserve(numPoints);//设置容器大小

  std::vector<int> numParentNeigh, numGrandParentNeigh;//创建容器存放父节点和祖父节点的个数
  numParentNeigh.resize(numPoints);//设置容器大小
  numGrandParentNeigh.resize(numPoints);

  // quant layer selection
  auto qpLayer = 0;

  // descend tree从根节点开始划分树
  weightsLf.resize(1);//将左侧树的大小设置为1
  attrsLf.resize(numAttrs);//将存放属性的容器大小设置为1
  for (int level = levelHfPos.size() - 1, isFirst = 1; level > 0; /*nop*/) {//获取当前层数信息
    int numNodes = weightsHf.size() - levelHfPos[level];//当前点云中所记录的所有重复点的个数减—当前层重复点的个数=当前层的重复点个数
    weightsLf.resize(weightsLf.size() + numNodes);//将左侧树的大小设置为当前层包含的点数
    attrsLf.resize(attrsLf.size() + numNodes * numAttrs);//将属性树的左侧设置为当前层包含的点数*属性值的个数
    //从顶层构建树
	expandLevel(
      level, numNodes, numAttrs, &weightsLf, &weightsHf, &attrsLf, &attrsHf);
    weightsHf.resize(levelHfPos[level]);//下一层中重复点的个数
    attrsHf.resize(levelHfPos[level] * numAttrs);//下一层中存放重复点的属性信息的容器的个数

    // expansion of level is complete, processing is now on the next level
    level--;//计算下一层

    // every three levels, perform transform
    //三次划分以后开始进行RAHT变换
    if (level % 3)//每三层进行一次RAHT变换，数学上表现为如果level可以整除3，则对当前左侧树中的点进行RAHT变换
      continue;

    // initial scan position of the coefficient buffer
    //  -> first level = all coeffs
    //  -> otherwise = ac coeffs only
    bool inheritDc = !isFirst;////判断是否继承父层DC系数
    bool enablePredictionInLvl = inheritDc && raht_prediction_enabled_flag;//属性预测开关，由inheritDc和raht_prediction_enabled_flag共同决定
    isFirst = 0;//其余变换打开属性预测开关

    // select quantiser according to transform layer
    qpLayer = std::min(qpLayer + 1, int(qpset.layers.size()) - 1);//选择量化参数

    // prepare reconstruction buffers
    //  previous reconstruction -> attrRecParent
	//将重建后各点的属性信息作为下一次变换中各点的父节点
    std::swap(attrRec, attrRecParent);
    std::swap(attrRecUs, attrRecParentUs);
    std::swap(numParentNeigh, numGrandParentNeigh);
    auto attrRecParentUsIt = attrRecParentUs.cbegin();
    auto attrRecParentIt = attrRecParent.cbegin();
    auto weightsParentIt = weightsParent.cbegin();
    auto numGrandParentNeighIt = numGrandParentNeigh.cbegin();

    for (int i = 0, iLast, iEnd = weightsLf.size(); i < iEnd; i = iLast) {
		//依次遍历每一个2*2*2的立方体
      // todo(df): hoist and dynamically allocate
      FixedPoint transformBuf[6][8] = {};//建立数组存放变换前后的系数
      FixedPoint(*transformPredBuf)[8] = &transformBuf[numAttrs];//创建指针读取变换系数
      int weights[8 + 8 + 8 + 8] = {};//存放权重信息
      Qps nodeQp[8] = {};//存放量化信息
      uint8_t occupancy = 0;//标记当前节点是否被占据

      // generate weights, occupancy mask, and fwd transform buffers
      // for all siblings of the current node.
	  //对第一个立方体（最多包含8个节点）进行预测和变换
      for (iLast = i; iLast < iEnd; iLast++) {
        int nextNode = iLast > i
          && !isSibling(weightsLf[iLast].pos, weightsLf[i].pos, level + 3);//判断当前点是否属于父块中的点
        if (nextNode)//如果不属于父节点，则跳出该次循环，判断下一个父节点和该节点的关系
          break;

		//若该节点属于当前父节点
        int nodeIdx = (weightsLf[iLast].pos >> level) & 0x7;//判断该节点在父节点中的位置
        weights[nodeIdx] = weightsLf[iLast].weight;//建立该父块中各子节点的权重索引
        nodeQp[nodeIdx][0] = weightsLf[iLast].qp[0] >> regionQpShift;//量化参数索引
        nodeQp[nodeIdx][1] = weightsLf[iLast].qp[1] >> regionQpShift;//量化参数索引

        occupancy |= 1 << nodeIdx;//计算当前父节点的占位码：用一个8ibt的二进制数标记父块中的八个子节点的占据情况

        if (isEncoder) {//编码
          for (int k = 0; k < numAttrs; k++)
            transformBuf[k][nodeIdx] = attrsLf[iLast * numAttrs + k];//将属性信息值存放在transformBuf中，同时左移15位，即保留15位小数点精度
        }
      }

      mkWeightTree(weights);//建立父节点内各个子节点的权重树：沿三个方向将相邻两个节点的权重进行相加自下而上得到变换后节点的权重大小

      if (!inheritDc) {//直接编码
        for (int j = i, nodeIdx = 0; nodeIdx < 8; nodeIdx++) {
          if (!weights[nodeIdx])//权重信息为0，判断下一个点
            continue;
          numParentNeigh[j++] = 19;//给它的每一个父邻居节点的个数赋初始值为19
        }
      }

      // Inter-level prediction:
      //  - Find the parent neighbours of the current node
      //  - Generate prediction for all attributes into transformPredBuf
      //  - Subtract transformed coefficients from forward transform
      //  - The transformPredBuf is then used for reconstruction
	  //属性信息预测：利用父节点和父节点的邻居信息，预测子节点的属性信息
      bool enablePrediction = enablePredictionInLvl;
      if (enablePredictionInLvl) {// 进行属性预测
        // indexes of the neighbouring parents
        int parentNeighIdx[19];//标记父邻居节点是否存在，-1表示不存在，标记邻居父节点相对于子节点的位置信息
        int parentNeighWeights[19];//邻居父节点的权重信息

        int parentNeighCount = 0;//初始邻居父节点个数为0
        if (*numGrandParentNeighIt < predictionThreshold[0]) {//邻居祖父节点的个数小于阈值，此处设置为2
          enablePrediction = false;//不进行预测
        } else {//寻找邻居父节点
          findNeighbours(
            weightsParent.cbegin(), weightsParent.cend(), weightsParentIt,
            level + 3, occupancy, parentNeighIdx, parentNeighWeights);
          for (int i = 0; i < 19; i++) {
            parentNeighCount += (parentNeighIdx[i] != -1);//统计找到的邻居数目：计算父节点的个数
          }
          if (parentNeighCount < predictionThreshold[1]) {//邻居父节点的个数小于阈值，不进行预测，此处阈值设置为4
            enablePrediction = false;//关闭预测开关
          } else//进行属性信息预测，计算预测值
          //利用父节点的邻居属性值对父节点中的每一个子节点属性值进行预测
		 intraDcPred(
              numAttrs, parentNeighIdx, parentNeighWeights, occupancy,
              attrRecParent.begin(), transformPredBuf);
        }

        for (int j = i, nodeIdx = 0; nodeIdx < 8; nodeIdx++) {
          if (!weights[nodeIdx])
            continue;
          numParentNeigh[j++] = parentNeighCount;//记录当前块的邻居节点个数
        }
      }

      int parentWeight = 0;//邻居父节点的权重信息
      if (inheritDc) {
        parentWeight = weightsParentIt->weight;//更新邻居父节点的权重信息
        weightsParentIt++;//指针指向下一个父块
        numGrandParentNeighIt++;//指向下一个父块
      }

      // normalise coefficients
	  //对属性信息值进行归一化，计算其属性均值
      for (int childIdx = 0; childIdx < 8; childIdx++) {
        if (weights[childIdx] <= 1)//判断当前子节点的权重小于等于1
          continue;//跳出该次循环，进行下一次判断

        // Summed attribute values
        if (isEncoder) {//编码
          FixedPoint rsqrtWeight;
          uint64_t w = weights[childIdx];
          int shift = (w > 1024 ? 5 : 0) + (w > 16384 ? 2 : 0)
            + (w > 262144 ? 2 : 0) + (w > 4194304 ? 2 : 0);
          rsqrtWeight.val = irsqrt(w) >> (40 - shift - FixedPoint::kFracBits);
          for (int k = 0; k < numAttrs; k++) {
            transformBuf[k][childIdx].val >>= shift;
            transformBuf[k][childIdx] *= rsqrtWeight;//属性信息值/权重开根号
          }
        }

        // Predicted attribute values
        if (enablePrediction) {//如果属性预测开关打开
          FixedPoint sqrtWeight;
          sqrtWeight.val =
            isqrt(uint64_t(weights[childIdx]) << (2 * FixedPoint::kFracBits));
          for (int k = 0; k < numAttrs; k++)
            transformPredBuf[k][childIdx] *= sqrtWeight;//属性预测值*权重开根号
        }
      }

      // forward transform:
      //  - encoder: transform both attribute sums and prediction
      //  - decoder: just transform prediction
      if (isEncoder && enablePrediction)//预测和编码同时进行，对属性值和属性预测值都进行RAHT变换
        fwdTransformBlock222<RahtKernel>(2 * numAttrs, transformBuf, weights);
      else if (isEncoder)//只编码不预测，对属性信息进行RAHT变换
        fwdTransformBlock222<RahtKernel>(numAttrs, transformBuf, weights);
      else if (enablePrediction)//不编码只预测，对属性预测信息进行RAHT变换
        fwdTransformBlock222<RahtKernel>(numAttrs, transformPredBuf, weights);

      // per-coefficient operations:
      //  - subtract transform domain prediction (encoder)
      //  - write out/read in quantised coefficients
      //  - inverse quantise + add transform domain prediction
      scanBlock(weights, [&](int idx) {
        // skip the DC coefficient unless at the root of the tree
		  //对DC系数不进行量化
        if (inheritDc && !idx)
          return;

        // subtract transformed prediction (skipping DC)
        if (isEncoder && enablePrediction) {
          for (int k = 0; k < numAttrs; k++) {
            transformBuf[k][idx] -= transformPredBuf[k][idx];//计算属性信息值和属性预测信息进行RAHT变换之后的残差信息
          }
        }

        // The RAHT transform
		//量化
		//不进行预测时，直接对属性信息进行量化
		//进行预测时，对预测残差进行量化
        auto quantizers = qpset.quantizers(qpLayer, nodeQp[idx]);
        for (int k = 0; k < numAttrs; k++) {
          auto& q = quantizers[std::min(k, int(quantizers.size()) - 1)];

          if (isEncoder) {//编码
            auto coeff = transformBuf[k][idx].round();
            assert(coeff <= INT_MAX && coeff >= INT_MIN);
            *coeffBufItK[k]++ = coeff =
              q.quantize(coeff << kFixedPointAttributeShift);//量化
            transformPredBuf[k][idx] +=
              divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift);//反量化
          } else {//解码
            int64_t coeff = *coeffBufItK[k]++;
            transformPredBuf[k][idx] +=
              divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift);//反量化
          }
        }
      });

	  //属性信息重建
      // replace DC coefficient with parent if inheritable
      if (inheritDc) {//继承DC系数
        for (int k = 0; k < numAttrs; k++) {//通过父节点属性信息得到当前块的DC系数
          attrRecParentIt++;
          int64_t val = *attrRecParentUsIt++;//父节点属性信息
          if (val > 0)
            transformPredBuf[k][0].val = val << (15 - 2);//赋给当前块的DC系数
          else
            transformPredBuf[k][0].val = -((-val) << (15 - 2));
        }
      }

      invTransformBlock222<RahtKernel>(numAttrs, transformPredBuf, weights);//RAHT反变换：对反量化系数进行逆变换得到重建属性值

      for (int j = i, nodeIdx = 0; nodeIdx < 8; nodeIdx++) {
        if (!weights[nodeIdx])//权重为0，则直接进行下一次判断
          continue;

        for (int k = 0; k < numAttrs; k++) {
          FixedPoint temp = transformPredBuf[k][nodeIdx];//读取当前的变换系数信息
          temp.val <<= 2;//右移两位
          attrRecUs[j * numAttrs + k] = temp.round();//左移15位存入属性重建信息
        }

        // scale values for next level 对重建属性值进行归一化
		//属性信息/根号weight，右移15位，舍去精度，得到当前父级块的平均属性信息
        if (weights[nodeIdx] > 1) {
          FixedPoint rsqrtWeight;
          uint64_t w = weights[nodeIdx];
          int shift = (w > 1024 ? 5 : 0) + (w > 16384 ? 2 : 0)
            + (w > 262144 ? 2 : 0) + (w > 4194304 ? 2 : 0);
          rsqrtWeight.val = irsqrt(w) >> (40 - shift - FixedPoint::kFracBits);
          for (int k = 0; k < numAttrs; k++) {
            transformPredBuf[k][nodeIdx].val >>= shift;
            transformPredBuf[k][nodeIdx] *= rsqrtWeight;//属性信息/weight开根号
          }
        }

        for (int k = 0; k < numAttrs; k++)
          attrRec[j * numAttrs + k] = transformPredBuf[k][nodeIdx].round();//将重建的平均属性信息值写入属性重建信息
        j++;
      }
    }

    // preserve current weights/positions for later search
    weightsParent = weightsLf;//该层中点的个数=下层中父节点的个数
  }

  // process duplicate points at level 0
  std::swap(attrRec, attrRecParent);
  auto attrRecParentIt = attrRecParent.cbegin();
  auto attrsHfIt = attrsHf.cbegin();

  for (int i = 0, out = 0, iEnd = weightsLf.size(); i < iEnd; i++) {
    int weight = weightsLf[i].weight;
    Qps nodeQp = {weightsLf[i].qp[0] >> regionQpShift,
                  weightsLf[i].qp[1] >> regionQpShift};

    // unique points have weight = 1
    if (weight == 1) {
      for (int k = 0; k < numAttrs; k++)
        attrRec[out++] = *attrRecParentIt++;
      continue;
    }

    // duplicates
    FixedPoint attrSum[3];
    FixedPoint attrRecDc[3];
    FixedPoint sqrtWeight;
    sqrtWeight.val = isqrt(uint64_t(weight) << (2 * FixedPoint::kFracBits));
    for (int k = 0; k < numAttrs; k++) {
      if (isEncoder)
        attrSum[k] = attrsLf[i * numAttrs + k];
      attrRecDc[k] = *attrRecParentIt++;
      attrRecDc[k] *= sqrtWeight;
    }

    FixedPoint rsqrtWeight;
    for (int w = weight - 1; w > 0; w--) {
      RahtKernel kernel(w, 1);
      int shift = (w > 1024 ? 5 : 0) + (w > 16384 ? 2 : 0)
        + (w > 262144 ? 2 : 0) + (w > 4194304 ? 2 : 0);
      if (isEncoder)
        rsqrtWeight.val = irsqrt(w) >> (40 - shift - FixedPoint::kFracBits);

      auto quantizers = qpset.quantizers(qpLayer, nodeQp);
      for (int k = 0; k < numAttrs; k++) {
        auto& q = quantizers[std::min(k, int(quantizers.size()) - 1)];

        FixedPoint transformBuf[2];
        if (isEncoder) {
          // invert the initial reduction (sum)
          // NB: read from (w-1) since left side came from attrsLf.
          transformBuf[1] = attrsHfIt[(w - 1) * numAttrs + k];
          attrSum[k] -= transformBuf[1];
          transformBuf[0] = attrSum[k];

          // NB: weight of transformBuf[1] is by construction 1.
          transformBuf[0].val >>= shift;
          transformBuf[0] *= rsqrtWeight;

          kernel.fwdTransform(
            transformBuf[0], transformBuf[1], &transformBuf[0],
            &transformBuf[1]);

          auto coeff = transformBuf[1].round();
          assert(coeff <= INT_MAX && coeff >= INT_MIN);
          *coeffBufItK[k]++ = coeff =
            q.quantize(coeff << kFixedPointAttributeShift);
          transformBuf[1] =
            divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift);
        } else {
          int64_t coeff = *coeffBufItK[k]++;
          transformBuf[1] =
            divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift);
        }

        // inherit the DC value
        transformBuf[0] = attrRecDc[k];

        kernel.invTransform(
          transformBuf[0], transformBuf[1], &transformBuf[0],
          &transformBuf[1]);

        attrRecDc[k] = transformBuf[0];
        attrRec[out + w * numAttrs + k] = transformBuf[1].round();
        if (w == 1)
          attrRec[out + k] = transformBuf[0].round();
      }
    }

    attrsHfIt += (weight - 1) * numAttrs;
    out += weight * numAttrs;
  }

  // write-back reconstructed attributes
  assert(attrRec.size() == numAttrs * numPoints);
  std::copy(attrRec.begin(), attrRec.end(), attributes);
}

//============================================================================
/*
 * RAHT Fixed Point
 *
 * Inputs:
 * quantStepSizeLuma = Quantization step
 * mortonCode = list of 'voxelCount' Morton codes of voxels, sorted in ascending Morton code order
 * attributes = 'voxelCount' x 'attribCount' array of attributes, in row-major order
 * attribCount = number of attributes (e.g., 3 if attributes are red, green, blue)
 * voxelCount = number of voxels
 *
 * Outputs:
 * weights = list of 'voxelCount' weights associated with each transform coefficient
 * coefficients = quantized transformed attributes array, in column-major order
 * binaryLayer = binary layer where each coefficient was generated
 *
 * Note output weights are typically used only for the purpose of
 * sorting or bucketing for entropy coding.
 */
void
regionAdaptiveHierarchicalTransform(
  bool raht_prediction_enabled_flag,
  const int predictionThreshold[2],
  const QpSet& qpset,
  const Qps* pointQpOffsets,
  int64_t* mortonCode,
  int* attributes,
  const int attribCount,
  const int voxelCount,
  int* coefficients)
{
  uraht_process<true>(
    raht_prediction_enabled_flag, predictionThreshold, qpset, pointQpOffsets,
    voxelCount, attribCount, mortonCode, attributes, coefficients);
}

//============================================================================
/*
 * inverse RAHT Fixed Point
 *
 * Inputs:
 * quantStepSizeLuma = Quantization step
 * mortonCode = list of 'voxelCount' Morton codes of voxels, sorted in ascending Morton code order
 * attribCount = number of attributes (e.g., 3 if attributes are red, green, blue)
 * voxelCount = number of voxels
 * coefficients = quantized transformed attributes array, in column-major order
 *
 * Outputs:
 * attributes = 'voxelCount' x 'attribCount' array of attributes, in row-major order
 *
 * Note output weights are typically used only for the purpose of
 * sorting or bucketing for entropy coding.
 */
void
regionAdaptiveHierarchicalInverseTransform(
  bool raht_prediction_enabled_flag,
  const int predictionThreshold[2],
  const QpSet& qpset,
  const Qps* pointQpOffsets,
  int64_t* mortonCode,
  int* attributes,
  const int attribCount,
  const int voxelCount,
  int* coefficients)
{
  uraht_process<false>(
    raht_prediction_enabled_flag, predictionThreshold, qpset, pointQpOffsets,
    voxelCount, attribCount, mortonCode, attributes, coefficients);
}

//============================================================================

}  // namespace pcc
