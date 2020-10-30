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
 //   numnodes:��ǰ�����������ĵ���
 //   numAttrs:ÿ���������������Ϣ�ĸ���������rgb��Ϊ3����Ϣ
 //   weightsIn:�������������еĵ�İ�����λ����Ϣ��Ȩ����Ϣ������������Ϣ
 //   weightsOut:������ǰ��ĵ����Ϣ
 //   attrsIn:�������������еĵ��������Ϣ
	//attrsOut:��ǰ���������Ϣ

  int64_t posPrev = -1; //����һ����ʼ����Ϊ-1
  auto weightsInWrIt = weightsIn->begin();//ָ����������еĵ���׵�ַ
  auto weightsInRdIt = weightsIn->cbegin();//ָ����������еĵ���׵�ַ��cbegin�����޸�ָ��ָ���Ԫ��
  auto attrsInWrIt = attrsIn->begin();//ָ�����������Ϣ���׵�ַ
  auto attrsInRdIt = attrsIn->begin();//ָ�����������Ϣ���׵�ַ
  for (int i = 0; i < numNodes; i++) {
    const auto& node = *weightsInRdIt++;//ͨ��node��¼��ǰ�ڵ����Ϣ

    // copy across unique nodes
    if (node.pos != posPrev) {//��ǰ��û���ظ����֣���Ī������Ϣ����һ���㲻ͬ
      posPrev = node.pos;//��¼��ǰ���Ī������Ϣ
      *weightsInWrIt++ = node;//�����ǰ�����Ϣ
      for (int k = 0; k < numAttrs; k++)
        *attrsInWrIt++ = *attrsInRdIt++;//�����ǰ���������Ϣ
      continue;
    }

    // duplicate node
	//��ǰ���ظ����֣�����ǰ������һ�����Ī������Ϣ��ͬ����ʱ�����ȥ���ظ���Ĳ���
    (weightsInWrIt - 1)->weight += node.weight;//������һ�����Ȩ����Ϣ=��ǰ���Ȩ����Ϣ+��һ�����Ȩ����Ϣ
    weightsOut->push_back(node);//��weightsOut�������ǰ�����Ϣ
    for (int k = 0; k < numAttrs; k++) {
      *(attrsInWrIt - numAttrs + k) += *attrsInRdIt;//������һ�����������Ϣֵ=��ǰ�������ֵ+��һ�����Ӧ������ֵ
      attrsOut->push_back(*attrsInRdIt++);//��attrsOut�������ǰ���������Ϣֵ
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
	// ͨ������ǰ�����һ���Ī����������ͬ��λ�����ж��������ڵ��Ƿ�����ͬһ���ڵ㣬�Ӷ����е�ĺϲ�
  //   level:��ǰ��Ĳ���
  //   numnodes:��һ�����������ĵ���
  //   numAttrs:ÿ���������������Ϣ�ĸ���������rgb��Ϊ3����Ϣ
  //   weightsIn:��һ���е�λ����Ϣ������������֮�����Ϊ��ǰ�������еĵ�İ�����λ����Ϣ��Ȩ����Ϣ������������Ϣ
  //   weightsOut:ͬһ���ڵ�ͬ��/�е���һ�������Ϣ
  //   attrsIn:����������һ����������Ϣ������Ϊ��ǰ��ĵ��������Ϣ
  //   attrsOut:ͬһ���ڵ�ͬ��/�е���һ�����������Ϣ
  int64_t posPrev = -1;//����һ����ʼ����Ϊ-1
  auto weightsInWrIt = weightsIn->begin();//ָ����������еĵ���׵�ַ
  auto weightsInRdIt = weightsIn->cbegin();//ָ����������еĵ���׵�ַ��cbegin�����޸�ָ��ָ���Ԫ��
  auto attrsInWrIt = attrsIn->begin();//ָ�����������Ϣ���׵�ַ
  auto attrsInRdIt = attrsIn->begin();//ָ�����������Ϣ���׵�ַ
  for (int i = 0; i < numNodes; i++) {//����Ī����˳����������еĵ�
    auto& node = *weightsInRdIt++;//node��¼��ǰ�����Ϣ
    bool newPair = (posPrev ^ node.pos) >> level != 0;
	//�жϵ�ǰ���Ī�������С��ǰһ�����Ƿ����ھӣ�newPairΪ��ʱ�����ھӣ�Ϊ��ʱ���ھ�
    posPrev = node.pos;
    posPrev = node.pos;//posPrev��¼��ǰ���Ī������Ϣ
    //��ǰ����ǰһ���㲻Ϊ�ھӹ�ϵʱ������ǰ�㼸����Ϣ����weightsLf��������Ϣ����attrsLf
    if (newPair) {
      *weightsInWrIt++ = node;//�����ǰ�����Ϣ��weightsIn��
      for (int k = 0; k < numAttrs; k++)
        *attrsInWrIt++ = *attrsInRdIt++;//�����ǰ���������Ϣ��attrsIn��
    } else {
      //���ھ�ʱ�����������ڵ��Ȩ����Ӵ���weightsLf��������Ӵ���attrsLf��Ī����Ϊǰһ�����Ī����ֵ������ǰ�㼸����Ϣ����weightsHf,������Ϣ����attrsHf
      auto& left = *(weightsInWrIt - 1);//��¼��ǰ�����һ���ڵ����Ϣ
      left.weight += node.weight;//������һ���ڵ��Ȩ����Ϣ=��һ���ڵ��Ȩ����Ϣ+��ǰ�ڵ��Ȩ����Ϣ
      left.qp[0] = (left.qp[0] + node.qp[0]) >> 1;//������һ���ڵ������������Ϣ=��ǰ��+��һ���ڵ��Ȩ����Ϣ/2
      left.qp[1] = (left.qp[1] + node.qp[1]) >> 1;//������һ���ڵ������������Ϣ=��ǰ��+��һ���ڵ��Ȩ����Ϣ/2
      weightsOut->push_back(node);//�����ǰ�㵽weightOut��

      for (int k = 0; k < numAttrs; k++) {
        *(attrsInWrIt - numAttrs + k) += *attrsInRdIt;//������һ�����������Ϣֵ=��ǰ�������ֵ+��һ�����Ӧ������ֵ
        attrsOut->push_back(*attrsInRdIt++);//�����ǰ�������ֵ��attrsOut��
      }
    }
  }

  // number of nodes in next level ���ظò�û���ھӵĵ���
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
  // ͨ������ǰ�����һ���Ī����������ͬ��λ�����ж��������ڵ��Ƿ�����ͬһ���ڵ㣬�Ӷ����е�ĺϲ�
  //   level:��ǰ��Ĳ���
  //   numnodes:��ǰ�������������ظ������ ��ǰ��ĵ���=��һ���е�ĸ���+��ǰ�����ظ���ĸ���
  //   numAttrs:ÿ���������������Ϣ�ĸ���������rgb��Ϊ3����Ϣ
  //   weightsIn:��ǰ�������еĵ�İ�����λ����Ϣ��Ȩ����Ϣ������������Ϣ
  //   weightsOut:��ǰ�����ظ������Ϣ
  //   attrsIn:��ǰ��ĵ��������Ϣ
  //   attrsOut:��ǰ�����ظ����������Ϣ
  if (numNodes == 0)//�жϸò��е����Ƿ�Ϊ0����Ϊ0����ֱ���˳��ú���
    return;

  // process a single level of the tree
  auto weightsInWrIt = weightsIn->rbegin();//��ȡ��ǰ�����һ�������Ϣ
  auto weightsInRdIt = std::next(weightsIn->crbegin(), numNodes);//��ȡ�ϲ����һ�������Ϣ
  auto weightsOutRdIt = weightsOut->crbegin();//�����ظ��������һ���ظ������Ϣ
  auto attrsInWrIt = attrsIn->rbegin();//��ǰ�����һ��������һ��������Ϣ
  auto attrsInRdIt = std::next(attrsIn->crbegin(), numNodes * numAttrs);//�������һ����ĵ�һ��������Ϣֵ
  auto attrsOutRdIt = attrsOut->crbegin();//�ظ��������һ���ظ�������һ��������Ϣֵ
  for (int i = 0; i < numNodes;) {
    //�жϵ�ǰ����weightsHf�дӺ�ʼ�����ĵ��Ƿ����ھӹ�ϵ��isPairΪ���ʾ���ھӹ�ϵ��Ϊ�ٱ�ʾ�����ھӹ�ϵ
    bool isPair = (weightsOutRdIt->pos ^ weightsInRdIt->pos) >> level == 0;
    if (!isPair) {//����������ھӽڵ�
      *weightsInWrIt++ = *weightsInRdIt++;//����ǰ������һ������Ϣ������ǰ�ظ���
      for (int k = 0; k < numAttrs; k++)
        *attrsInWrIt++ = *attrsInRdIt++;//����ǰ������һ����������Ϣ����ǰ���ظ���
      continue;
    }

	//�����ھӽڵ㣬��ǰ������Լ�ȥattrsHf�Ӻ�ʼ�����������ֵ�õ���ǰ���ھӽڵ������ֵ
    // going to process a pair
    i++;//�ظ������+1

    // Out node is inserted before In node.
    const auto& nodeDelta = *weightsInWrIt++ = *weightsOutRdIt++;//����ǰ�ظ������Ϣ������ǰ��ĵ�
    auto curAttrIt = attrsInWrIt;//��ȡ��ǰ���е������ֵ
    for (int k = 0; k < numAttrs; k++)
      *attrsInWrIt++ = *attrsOutRdIt++;//�ظ����ֵ������ǰ��

    // move In node to correct position, subtracting delta
    *weightsInWrIt = *weightsInRdIt++;//��һ���ж�Ӧ���λ��
    (weightsInWrIt++)->weight -= nodeDelta.weight;//�����һ���е��ڵ�ǰ���Ȩ��
    for (int k = 0; k < numAttrs; k++) {
      *attrsInWrIt = *attrsInRdIt++;//��һ���ж�Ӧ���������Ϣ
      *attrsInWrIt++ -= *curAttrIt++;//�����һ���е��ڵ�ǰ���������Ϣ
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
	//1-19�ֱ��ǹ̶�λ�õ��ھӽڵ㣬[1]��ʾ���������棬���ߵ��ھ�
  //���ڵ�ǰ�ڵ㣬���ھ���6�����棬12�����߰�����ǰ���ڵ㱾��һ�������19���ھӣ�neighMasks��ʾ��ǰ���ڵ���
  //���ھӹ�����ߵ���Ч�ӽڵ㣬1��ʾ��Ч��0��ʾ��Ч������15��00001111������ʾ���������ھӣ���ǰ���ڵ�������Ϊ0,1,2,3���ӽڵ�����Ч�ӽڵ㡣
  static const uint8_t neighMasks[19] = {255, 15, 240, 51, 204, 85,  170,
                                         3,   12, 5,   10, 48,  192, 80,
                                         160, 17, 34,  68, 136};
  //��8bit��ʾ��ǰ�ڵ��Ԥ����ӽڵ��λ�ã�8bit�ֱ��ʾ��ǰ�������ӽڵ��λ����Ϣ��1��ʾ�ɽ���Ԥ�⣬0��ʾ����Ԥ����ӿ�

  // current position (discard extra precision)
  int64_t cur_pos = it->pos >> level;//���㵱ǰ���ڵ�ĵ�ַ

  // the position of the parent, offset by (-1,-1,-1)
  int64_t base_pos = morton3dAdd(cur_pos, -1ll);//�Ե�ǰ���ڵ��������ƫ��

  // these neighbour offsets are relative to base_pos
  static const uint8_t neighOffset[19] = {0,  3,  35, 5,  21, 6, 14, 1,  17, 2,
                                          10, 33, 49, 34, 42, 4, 12, 20, 28};

   //��ʾÿһ���ھӽڵ�����ڵ�ǰ�ڵ�ļ�������ƫ����
  // special case for the direct parent (no need to search);
  parentNeighIdx[0] = std::distance(first, it);
  parentNeighWeights[0] = it->weight;//��0���ھ�Ϊ��ǰ���ڵ㱾��

  for (int i = 1; i < 19; i++) {
    // Only look for neighbours that have an effect
    if (!(occupancy & neighMasks[i])) {
      parentNeighIdx[i] = -1;
      continue;
    }
    // �������ڵ��ڶ��ھ���Ч���ӽڵ�ռ��ʱ�Ż��������ھӣ����統����Ϊ0, 1, 2,3���ӽڵ�Ϊ��ʱ������������๲���ھ�

    // compute neighbour address to look for
    // the delta between it and the current position is
    int64_t neigh_pos = morton3dAdd(base_pos, neighOffset[i]);//�����ھӽڵ�Ը��ڵ������ƫ���������ھӽڵ��Ī����ֵ
    int64_t delta = neigh_pos - cur_pos;//�����ھӽڵ�����ڵ�ǰ���ڵ��λ�òв�


    // find neighbour
    auto found = findNeighbour(
      first, last, it, neigh_pos, delta,
      [=](decltype(*it)& candidate, int64_t neigh_pos) {
        return (candidate.pos >> level) < neigh_pos;
      });//�����ھ��뵱ǰ���ڵ�ľ���в������ھ�

    if (found == last) {
      parentNeighIdx[i] = -1;
      continue;
    }//���ҵ����ھӵ�ַ������ʱ�����ھӲ����ڣ�������Ϊ-1

    if ((found->pos >> level) != neigh_pos) {
      parentNeighIdx[i] = -1;
      continue;
    }//�����ݾ���в��ҵ����ھ�Ī��������level�Ժ󲻵����ھӵ�ַ�����ھӲ����ڣ�������Ϊ-1

    parentNeighIdx[i] = std::distance(first, found);//�����ھ�����
    parentNeighWeights[i] = found->weight;//�����ھ�Ȩ��
  }//��������18���ھ�
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
  //ÿһ���ھ�ֻ�Ե�ǰ���ڵ������乲����ߵ��ӽڵ��������Ԥ�⣬��predMasks��ֵתΪ�������Ժ�Ϊ1���ӽڵ���ø��ھӽ���Ԥ�⣬Ϊ0���ӽڵ㲻����Ԥ�⣬
  //���磬predMasks[1]=15����Ϊ00001111����ʾ��๲���ھ�ֻ����Ԥ������Ϊ0,1,2,3���ӽڵ㣬��predMasks[0]=255��ʾ��ǰ���ڵ�Ԥ�������ӽڵ㣻
  static const uint8_t predMasks[19] = {255, 15, 240, 51, 204, 85,  170,
                                        3,   12, 5,   10, 48,  192, 80,
                                        160, 17, 34,  68, 136};
  //Ԥ��Ȩ��
  //Ԥ��Ȩ�أ���ǰ���ڵ���ӽڵ��Ԥ��Ȩ��Ϊ4�������ھ�Ԥ��Ȩ��Ϊ2�������ھ�Ԥ��Ȩ��Ϊ1
  static const int predWeight[19] = {4, 2, 2, 2, 2, 2, 2, 1, 1, 1,
                                     1, 1, 1, 1, 1, 1, 1, 1, 1};

    //�滻���������Ĳ�ѯ��
  static const int kDivisors[25] = {8192, 6554, 5461, 4681, 4096, 3641, 3277,
                                    2979, 2731, 2521, 2341, 2185, 2048, 1928,
                                    1820, 1725, 1638, 1560, 1489, 1425, 1365,
                                    1311, 1260, 1214, 1170};

  int weightSum[8] = {-4, -4, -4, -4, -4, -4, -4, -4};//Ԥ��Ȩ��֮��

  std::fill_n(&predBuf[0][0], 8 * numAttrs, FixedPoint(0));

  int64_t neighValue[3];
  //ȥ����Ч�ھӵ���ֵ
  int64_t limitLow = 0;//Ԥ����ֵ����
  int64_t limitHigh = 0;//Ԥ����ֵ����
  for (int i = 0; i < 19; i++) {
    if (neighIdx[i] == -1)//�õ㲻���ڸ��ھӽڵ�
      continue;//�˳��ô�ѭ��

    auto neighValueIt = std::next(first, numAttrs * neighIdx[i]);
    for (int k = 0; k < numAttrs; k++)
      neighValue[k] = *neighValueIt++;

    // skip neighbours that are outside of threshold limits
	//�ж��ھӸ��ڵ��Ƿ�����ֵ�ڣ������ڣ���õ㲻����ΪԤ��ֵȥԤ���ӿ����Ϣ
    //���ھ��뵱ǰ���ڵ������ֵ����С��limitLow�����limitHighʱ��ȥ�����ھӣ�����ֻ����Y����������ֵ
    if (i) {
      if (10 * neighValue[0] <= limitLow || 10 * neighValue[0] >= limitHigh)
        continue;
    } else {
      constexpr int ratioThreshold1 = 2;
      constexpr int ratioThreshold2 = 25;
      limitLow = ratioThreshold1 * neighValue[0];//�����ھ����Ե���С��ֵ
      limitHigh = ratioThreshold2 * neighValue[0];//�����ھ����Ե������ֵ
    }

    // apply weighted neighbour value to masked positions
    for (int k = 0; k < numAttrs; k++)
      neighValue[k] *= predWeight[i] << pcc::FixedPoint::kFracBits;

    //��Ȩ��Ԥ������ֵ
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
  //��Ԥ������ֵ���й�һ��
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
//RAHT�任
template<class Kernel>
void
fwdTransformBlock222(
  int numBufs, FixedPoint buf[][8], int weights[8 + 8 + 8 + 8])
{
  static const int a[4 + 4 + 4] = {0, 2, 4, 6, 0, 4, 1, 5, 0, 1, 2, 3};//RAHT�任˳��
  static const int b[4 + 4 + 4] = {1, 3, 5, 7, 2, 6, 3, 7, 4, 5, 6, 7};
  for (int i = 0, iw = 0; i < 12; i++, iw += 2) {
    int i0 = a[i];//������б任�Ŀ������
    int i1 = b[i];

    if (weights[iw] + weights[iw + 1] == 0)//�������ڵ��Ȩ��֮��Ϊ0����õ㲻����RAHT�任
      continue;

    // only one occupied, propagate to next level
	//�������RAHT�任�������ڵ㣬ֻ��һ���ڵ㱻ռ�ݣ���ôֱ�ӽ��������ڵ��λ�ã�������һ�㣬׼�������´α任
    if (!weights[iw] || !weights[iw + 1]) {
      if (!weights[iw]) {
        for (int k = 0; k < numBufs; k++)
          std::swap(buf[k][i0], buf[k][i1]);
      }
      continue;
    }

	//�����ڵ㶼��ռ�ݣ�����RAHT�任
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
  //raht_prediction_enabled_flag:����Ԥ�⿪��
  //predictionThreshold[2]:����Ԥ����ֵ���жϸ��ڵ�������游�ڵ����
  //qpset:��������
  //pointQpOffsets:����ƫ�Ʋ���
  //numPoints:���������еĵ���
  //numAttrs:������Ϣ����
  //positions:λ�ã���Ī������Ϣ
  //attributes:������Ϣ
  //coeffBufIt:�任ϵ����Ϣ
  // coefficients are stored in three planar arrays.  coeffBufItK is a set
  // of iterators to each array.
  //������������Ϣ��Ӧ�ı任ϵ����������׵�ַд��coeffBufItK[3]��
  int32_t* coeffBufItK[3] = {
    coeffBufIt,
    coeffBufIt + numPoints,
    coeffBufIt + numPoints * 2,
  };

  
  if (numPoints == 1) {//����ǰ���������еĵ���Ϊ1
    auto quantizers = qpset.quantizers(0, pointQpOffsets[0]);//ѡ������������Y������Cb��Cr֮����ܴ�������ƫ��
    for (int k = 0; k < numAttrs; k++) {
      auto& q = quantizers[std::min(k, int(quantizers.size()) - 1)];//ѡ������������Y����ѡ��quantizers[0],Cb��Cr����ѡ��quantizers[1];

      if (isEncoder) {//����
        auto coeff = attributes[k];//�ѵ�ǰ���������Ϣֵ����coeff
        assert(coeff <= INT_MAX && coeff >= INT_MIN);//�жϵ�ǰ������Ϣֵ��[-2^32,2^32]֮��
        *coeffBufItK[k]++ = coeff =
          q.quantize(coeff << kFixedPointAttributeShift);//�Ե�ǰ������Ϣֵcoeff����������������coeffBufItk��
        attributes[k] =
          divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift);//������֮�������ֵ���з�������д��attributes��
      } else {//����
        int64_t coeff = *coeffBufItK[k]++;//��ȡ����֮�����Ϣ
        attributes[k] =
          divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift);//���з�����
      }
    }
    return;
  }

  std::vector<UrahtNode> weightsLf, weightsHf;
  //����RAHT�任����������֦����������ΪUrahtNode���ͣ�������λ����Ϣ��Ī���룬Ȩ����Ϣ������������Ϣ
  //�ڹ���RAHT�任��ʱ�����������ϸ���״̬��Ϣ���Ҳ������¼�ظ������Ϣ
  std::vector<int> attrsLf, attrsHf;//�������������Ϣ��RAHT�任��������֦

  weightsLf.reserve(numPoints);//������֦�Ĵ�С=�������еĵ���
  attrsLf.reserve(numPoints * numAttrs);//��֦������֦�Ĵ�С=�������еĵ���*ÿ�������������ֵ�ĸ���

  int regionQpShift = 4;//����ԭʼ������ƫ����Ϊ4

  // copy positions into internal form
  // todo(df): lift to api
  for (int i = 0; i < numPoints; i++) {//�������֦�а�Ī������Ī������Ϣ��Ȩ����Ϣ������ƫ������Ϣ��ÿ����ĳ�ʼȨ������Ϊ1
    weightsLf.emplace_back(UrahtNode{positions[i],
                                     1,
                                     {pointQpOffsets[i][0] << regionQpShift,
                                      pointQpOffsets[i][1] << regionQpShift}});
    for (int k = 0; k < numAttrs; k++) {
      attrsLf.push_back(attributes[i * numAttrs + k]);//��RAHT������Ϣ�������֦���δ��ÿ����ĸ���������Ϣ
    }
  }

  weightsHf.reserve(numPoints);//������֦�Ĵ�С=�������еĵ���
  attrsHf.reserve(numPoints * numAttrs);//��֦������֦�Ĵ�С=�������еĵ���*ÿ�������������ֵ�ĸ���

  // ascend tree
  std::vector<int> levelHfPos;//����levelHfPos���ÿ�����а����ĵ�ĸ���
  //����RAHTԤ����

  for (int level = 0, numNodes = weightsLf.size(); numNodes > 1; level++) {
    levelHfPos.push_back(weightsHf.size());
    //�洢ÿһ����ǰһ�������ھӵĳɶԵĵ���
    if (level == 0) {//����ǰ��Ϊ��0��
      // process any duplicate points�ڵ�һ��ȥ��ԭʼ�����е��ظ���
	//ʹ�ú���reduceUnique������0����������ظò��а����ĵ�ĸ���
      numNodes = reduceUnique(
        numNodes, numAttrs, &weightsLf, &weightsHf, &attrsLf, &attrsHf);
    } else {
      // normal level reduction �ж�ÿһ�����ھӹ�ϵ
		//��㹹������RAHT�任����
      numNodes = reduceLevel(
        level, numNodes, numAttrs, &weightsLf, &weightsHf, &attrsLf, &attrsHf);
    }
  }

  assert(weightsLf[0].weight == numPoints);//��֤���������еĵ���󶼻�ϲ��ڵ�һ�����У�������Ȩ�ص������е��Ȩ��֮��

  // reconstruction buffers
  std::vector<int> attrRec, attrRecParent;//������������ؽ���Ϣ
  attrRec.resize(numPoints * numAttrs);//����������С
  attrRecParent.resize(numPoints * numAttrs);

  std::vector<int> attrRecUs, attrRecParentUs;//������������ؽ���Ϣ
  attrRecUs.resize(numPoints * numAttrs);//����������С
  attrRecParentUs.resize(numPoints * numAttrs);

  std::vector<UrahtNode> weightsParent;//����������Ÿ��ڵ��Ȩ����Ϣ
  weightsParent.reserve(numPoints);//����������С

  std::vector<int> numParentNeigh, numGrandParentNeigh;//����������Ÿ��ڵ���游�ڵ�ĸ���
  numParentNeigh.resize(numPoints);//����������С
  numGrandParentNeigh.resize(numPoints);

  // quant layer selection
  auto qpLayer = 0;

  // descend tree�Ӹ��ڵ㿪ʼ������
  weightsLf.resize(1);//��������Ĵ�С����Ϊ1
  attrsLf.resize(numAttrs);//��������Ե�������С����Ϊ1
  for (int level = levelHfPos.size() - 1, isFirst = 1; level > 0; /*nop*/) {//��ȡ��ǰ������Ϣ
    int numNodes = weightsHf.size() - levelHfPos[level];//��ǰ����������¼�������ظ���ĸ���������ǰ���ظ���ĸ���=��ǰ����ظ������
    weightsLf.resize(weightsLf.size() + numNodes);//��������Ĵ�С����Ϊ��ǰ������ĵ���
    attrsLf.resize(attrsLf.size() + numNodes * numAttrs);//�����������������Ϊ��ǰ������ĵ���*����ֵ�ĸ���
    //�Ӷ��㹹����
	expandLevel(
      level, numNodes, numAttrs, &weightsLf, &weightsHf, &attrsLf, &attrsHf);
    weightsHf.resize(levelHfPos[level]);//��һ�����ظ���ĸ���
    attrsHf.resize(levelHfPos[level] * numAttrs);//��һ���д���ظ����������Ϣ�������ĸ���

    // expansion of level is complete, processing is now on the next level
    level--;//������һ��

    // every three levels, perform transform
    //���λ����Ժ�ʼ����RAHT�任
    if (level % 3)//ÿ�������һ��RAHT�任����ѧ�ϱ���Ϊ���level��������3����Ե�ǰ������еĵ����RAHT�任
      continue;

    // initial scan position of the coefficient buffer
    //  -> first level = all coeffs
    //  -> otherwise = ac coeffs only
    bool inheritDc = !isFirst;////�ж��Ƿ�̳и���DCϵ��
    bool enablePredictionInLvl = inheritDc && raht_prediction_enabled_flag;//����Ԥ�⿪�أ���inheritDc��raht_prediction_enabled_flag��ͬ����
    isFirst = 0;//����任������Ԥ�⿪��

    // select quantiser according to transform layer
    qpLayer = std::min(qpLayer + 1, int(qpset.layers.size()) - 1);//ѡ����������

    // prepare reconstruction buffers
    //  previous reconstruction -> attrRecParent
	//���ؽ�������������Ϣ��Ϊ��һ�α任�и���ĸ��ڵ�
    std::swap(attrRec, attrRecParent);
    std::swap(attrRecUs, attrRecParentUs);
    std::swap(numParentNeigh, numGrandParentNeigh);
    auto attrRecParentUsIt = attrRecParentUs.cbegin();
    auto attrRecParentIt = attrRecParent.cbegin();
    auto weightsParentIt = weightsParent.cbegin();
    auto numGrandParentNeighIt = numGrandParentNeigh.cbegin();

    for (int i = 0, iLast, iEnd = weightsLf.size(); i < iEnd; i = iLast) {
		//���α���ÿһ��2*2*2��������
      // todo(df): hoist and dynamically allocate
      FixedPoint transformBuf[6][8] = {};//���������ű任ǰ���ϵ��
      FixedPoint(*transformPredBuf)[8] = &transformBuf[numAttrs];//����ָ���ȡ�任ϵ��
      int weights[8 + 8 + 8 + 8] = {};//���Ȩ����Ϣ
      Qps nodeQp[8] = {};//���������Ϣ
      uint8_t occupancy = 0;//��ǵ�ǰ�ڵ��Ƿ�ռ��

      // generate weights, occupancy mask, and fwd transform buffers
      // for all siblings of the current node.
	  //�Ե�һ�������壨������8���ڵ㣩����Ԥ��ͱ任
      for (iLast = i; iLast < iEnd; iLast++) {
        int nextNode = iLast > i
          && !isSibling(weightsLf[iLast].pos, weightsLf[i].pos, level + 3);//�жϵ�ǰ���Ƿ����ڸ����еĵ�
        if (nextNode)//��������ڸ��ڵ㣬�������ô�ѭ�����ж���һ�����ڵ�͸ýڵ�Ĺ�ϵ
          break;

		//���ýڵ����ڵ�ǰ���ڵ�
        int nodeIdx = (weightsLf[iLast].pos >> level) & 0x7;//�жϸýڵ��ڸ��ڵ��е�λ��
        weights[nodeIdx] = weightsLf[iLast].weight;//�����ø����и��ӽڵ��Ȩ������
        nodeQp[nodeIdx][0] = weightsLf[iLast].qp[0] >> regionQpShift;//������������
        nodeQp[nodeIdx][1] = weightsLf[iLast].qp[1] >> regionQpShift;//������������

        occupancy |= 1 << nodeIdx;//���㵱ǰ���ڵ��ռλ�룺��һ��8ibt�Ķ���������Ǹ����еİ˸��ӽڵ��ռ�����

        if (isEncoder) {//����
          for (int k = 0; k < numAttrs; k++)
            transformBuf[k][nodeIdx] = attrsLf[iLast * numAttrs + k];//��������Ϣֵ�����transformBuf�У�ͬʱ����15λ��������15λС���㾫��
        }
      }

      mkWeightTree(weights);//�������ڵ��ڸ����ӽڵ��Ȩ�������������������������ڵ��Ȩ�ؽ���������¶��ϵõ��任��ڵ��Ȩ�ش�С

      if (!inheritDc) {//ֱ�ӱ���
        for (int j = i, nodeIdx = 0; nodeIdx < 8; nodeIdx++) {
          if (!weights[nodeIdx])//Ȩ����ϢΪ0���ж���һ����
            continue;
          numParentNeigh[j++] = 19;//������ÿһ�����ھӽڵ�ĸ�������ʼֵΪ19
        }
      }

      // Inter-level prediction:
      //  - Find the parent neighbours of the current node
      //  - Generate prediction for all attributes into transformPredBuf
      //  - Subtract transformed coefficients from forward transform
      //  - The transformPredBuf is then used for reconstruction
	  //������ϢԤ�⣺���ø��ڵ�͸��ڵ���ھ���Ϣ��Ԥ���ӽڵ��������Ϣ
      bool enablePrediction = enablePredictionInLvl;
      if (enablePredictionInLvl) {// ��������Ԥ��
        // indexes of the neighbouring parents
        int parentNeighIdx[19];//��Ǹ��ھӽڵ��Ƿ���ڣ�-1��ʾ�����ڣ�����ھӸ��ڵ�������ӽڵ��λ����Ϣ
        int parentNeighWeights[19];//�ھӸ��ڵ��Ȩ����Ϣ

        int parentNeighCount = 0;//��ʼ�ھӸ��ڵ����Ϊ0
        if (*numGrandParentNeighIt < predictionThreshold[0]) {//�ھ��游�ڵ�ĸ���С����ֵ���˴�����Ϊ2
          enablePrediction = false;//������Ԥ��
        } else {//Ѱ���ھӸ��ڵ�
          findNeighbours(
            weightsParent.cbegin(), weightsParent.cend(), weightsParentIt,
            level + 3, occupancy, parentNeighIdx, parentNeighWeights);
          for (int i = 0; i < 19; i++) {
            parentNeighCount += (parentNeighIdx[i] != -1);//ͳ���ҵ����ھ���Ŀ�����㸸�ڵ�ĸ���
          }
          if (parentNeighCount < predictionThreshold[1]) {//�ھӸ��ڵ�ĸ���С����ֵ��������Ԥ�⣬�˴���ֵ����Ϊ4
            enablePrediction = false;//�ر�Ԥ�⿪��
          } else//����������ϢԤ�⣬����Ԥ��ֵ
          //���ø��ڵ���ھ�����ֵ�Ը��ڵ��е�ÿһ���ӽڵ�����ֵ����Ԥ��
		 intraDcPred(
              numAttrs, parentNeighIdx, parentNeighWeights, occupancy,
              attrRecParent.begin(), transformPredBuf);
        }

        for (int j = i, nodeIdx = 0; nodeIdx < 8; nodeIdx++) {
          if (!weights[nodeIdx])
            continue;
          numParentNeigh[j++] = parentNeighCount;//��¼��ǰ����ھӽڵ����
        }
      }

      int parentWeight = 0;//�ھӸ��ڵ��Ȩ����Ϣ
      if (inheritDc) {
        parentWeight = weightsParentIt->weight;//�����ھӸ��ڵ��Ȩ����Ϣ
        weightsParentIt++;//ָ��ָ����һ������
        numGrandParentNeighIt++;//ָ����һ������
      }

      // normalise coefficients
	  //��������Ϣֵ���й�һ�������������Ծ�ֵ
      for (int childIdx = 0; childIdx < 8; childIdx++) {
        if (weights[childIdx] <= 1)//�жϵ�ǰ�ӽڵ��Ȩ��С�ڵ���1
          continue;//�����ô�ѭ����������һ���ж�

        // Summed attribute values
        if (isEncoder) {//����
          FixedPoint rsqrtWeight;
          uint64_t w = weights[childIdx];
          int shift = (w > 1024 ? 5 : 0) + (w > 16384 ? 2 : 0)
            + (w > 262144 ? 2 : 0) + (w > 4194304 ? 2 : 0);
          rsqrtWeight.val = irsqrt(w) >> (40 - shift - FixedPoint::kFracBits);
          for (int k = 0; k < numAttrs; k++) {
            transformBuf[k][childIdx].val >>= shift;
            transformBuf[k][childIdx] *= rsqrtWeight;//������Ϣֵ/Ȩ�ؿ�����
          }
        }

        // Predicted attribute values
        if (enablePrediction) {//�������Ԥ�⿪�ش�
          FixedPoint sqrtWeight;
          sqrtWeight.val =
            isqrt(uint64_t(weights[childIdx]) << (2 * FixedPoint::kFracBits));
          for (int k = 0; k < numAttrs; k++)
            transformPredBuf[k][childIdx] *= sqrtWeight;//����Ԥ��ֵ*Ȩ�ؿ�����
        }
      }

      // forward transform:
      //  - encoder: transform both attribute sums and prediction
      //  - decoder: just transform prediction
      if (isEncoder && enablePrediction)//Ԥ��ͱ���ͬʱ���У�������ֵ������Ԥ��ֵ������RAHT�任
        fwdTransformBlock222<RahtKernel>(2 * numAttrs, transformBuf, weights);
      else if (isEncoder)//ֻ���벻Ԥ�⣬��������Ϣ����RAHT�任
        fwdTransformBlock222<RahtKernel>(numAttrs, transformBuf, weights);
      else if (enablePrediction)//������ֻԤ�⣬������Ԥ����Ϣ����RAHT�任
        fwdTransformBlock222<RahtKernel>(numAttrs, transformPredBuf, weights);

      // per-coefficient operations:
      //  - subtract transform domain prediction (encoder)
      //  - write out/read in quantised coefficients
      //  - inverse quantise + add transform domain prediction
      scanBlock(weights, [&](int idx) {
        // skip the DC coefficient unless at the root of the tree
		  //��DCϵ������������
        if (inheritDc && !idx)
          return;

        // subtract transformed prediction (skipping DC)
        if (isEncoder && enablePrediction) {
          for (int k = 0; k < numAttrs; k++) {
            transformBuf[k][idx] -= transformPredBuf[k][idx];//����������Ϣֵ������Ԥ����Ϣ����RAHT�任֮��Ĳв���Ϣ
          }
        }

        // The RAHT transform
		//����
		//������Ԥ��ʱ��ֱ�Ӷ�������Ϣ��������
		//����Ԥ��ʱ����Ԥ��в��������
        auto quantizers = qpset.quantizers(qpLayer, nodeQp[idx]);
        for (int k = 0; k < numAttrs; k++) {
          auto& q = quantizers[std::min(k, int(quantizers.size()) - 1)];

          if (isEncoder) {//����
            auto coeff = transformBuf[k][idx].round();
            assert(coeff <= INT_MAX && coeff >= INT_MIN);
            *coeffBufItK[k]++ = coeff =
              q.quantize(coeff << kFixedPointAttributeShift);//����
            transformPredBuf[k][idx] +=
              divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift);//������
          } else {//����
            int64_t coeff = *coeffBufItK[k]++;
            transformPredBuf[k][idx] +=
              divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift);//������
          }
        }
      });

	  //������Ϣ�ؽ�
      // replace DC coefficient with parent if inheritable
      if (inheritDc) {//�̳�DCϵ��
        for (int k = 0; k < numAttrs; k++) {//ͨ�����ڵ�������Ϣ�õ���ǰ���DCϵ��
          attrRecParentIt++;
          int64_t val = *attrRecParentUsIt++;//���ڵ�������Ϣ
          if (val > 0)
            transformPredBuf[k][0].val = val << (15 - 2);//������ǰ���DCϵ��
          else
            transformPredBuf[k][0].val = -((-val) << (15 - 2));
        }
      }

      invTransformBlock222<RahtKernel>(numAttrs, transformPredBuf, weights);//RAHT���任���Է�����ϵ��������任�õ��ؽ�����ֵ

      for (int j = i, nodeIdx = 0; nodeIdx < 8; nodeIdx++) {
        if (!weights[nodeIdx])//Ȩ��Ϊ0����ֱ�ӽ�����һ���ж�
          continue;

        for (int k = 0; k < numAttrs; k++) {
          FixedPoint temp = transformPredBuf[k][nodeIdx];//��ȡ��ǰ�ı任ϵ����Ϣ
          temp.val <<= 2;//������λ
          attrRecUs[j * numAttrs + k] = temp.round();//����15λ���������ؽ���Ϣ
        }

        // scale values for next level ���ؽ�����ֵ���й�һ��
		//������Ϣ/����weight������15λ����ȥ���ȣ��õ���ǰ�������ƽ��������Ϣ
        if (weights[nodeIdx] > 1) {
          FixedPoint rsqrtWeight;
          uint64_t w = weights[nodeIdx];
          int shift = (w > 1024 ? 5 : 0) + (w > 16384 ? 2 : 0)
            + (w > 262144 ? 2 : 0) + (w > 4194304 ? 2 : 0);
          rsqrtWeight.val = irsqrt(w) >> (40 - shift - FixedPoint::kFracBits);
          for (int k = 0; k < numAttrs; k++) {
            transformPredBuf[k][nodeIdx].val >>= shift;
            transformPredBuf[k][nodeIdx] *= rsqrtWeight;//������Ϣ/weight������
          }
        }

        for (int k = 0; k < numAttrs; k++)
          attrRec[j * numAttrs + k] = transformPredBuf[k][nodeIdx].round();//���ؽ���ƽ��������Ϣֵд�������ؽ���Ϣ
        j++;
      }
    }

    // preserve current weights/positions for later search
    weightsParent = weightsLf;//�ò��е�ĸ���=�²��и��ڵ�ĸ���
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
