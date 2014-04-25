/*!
  @file rx_ssm_tables.h
	
  @brief Marching Square用のテーブル

  @author Makoto Fujisawa
  @date 2011-05
*/
// FILE --rx_ssm_tables.h--

#ifndef _RX_SSM_TABLES_H_
#define _RX_SSM_TABLES_H_



// 三角形メッシュテーブル
//  - ノード, エッジ上の頂点の位置を示す8ビットフラグから
//    三角形のグリッド内頂点インデックスを示すテーブル
// 
// ノード頂点
// 3 - 2
// |   |
// 0 - 1
// 	 
// エッジ頂点
// - 6 -
// 7   5
// - 4 -
// 
// エッジ頂点(back vertex)
//  - 10 -
// 11     9
//  -  8 -
// 
// エッジ頂点2(back-2 vertex)
// 12, 13
//
#define X -1
// 上位4ビット : エッジ頂点
// 下位4ビット : ノード頂点
// * : 外部輪郭, + : 内部輪郭, *+ : 内部/外部輪郭混合
const int g_MeshTable[256][19] = {
	// パターン0 : 内部メッシュ
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0000 0000
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0000 0001
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0000 0010
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0000 0011
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0000 0100
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0000 0101
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0000 0110
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0000 0111
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0000 1000
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0000 1001
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0000 1010
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0000 1011
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0000 1100
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0000 1101
	{2, 0, 1, 2, 0, 2, 3, X, X, X, X, X, X, X, X, X, X, X, X},	// 0000 1110
	{2, 0, 1, 3, 3, 1, 2, X, X, X, X, X, X, X, X, X, X, X, X},	// 0000 1111 *

	// パターン1 : 内部輪郭(輪郭始点)
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0001 0000
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0001 0001
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0001 0010
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0001 0011
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0001 0100
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0001 0101
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0001 0110
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0001 0111
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0001 1000
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0001 1001
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0001 1010
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0001 1011
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0001 1100
	{3, 0, 4, 3, 3, 4, 2, 1, 2, 8, X, X, X, X, X, X, X, X, X},	// 0001 1101 +
	{3, 0, 8, 3, 3, 4, 2, 1, 2, 4, X, X, X, X, X, X, X, X, X},	// 0001 1110 +
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0001 1111

	// パターン2 : 内部輪郭(輪郭始点)
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0010 0000
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0010 0001
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0010 0010
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0010 0011
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0010 0100
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0010 0101
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0010 0110
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0010 0111
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0010 1000
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0010 1001
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0010 1010
	{3, 0, 1, 5, 0, 5, 3, 3, 9, 2, X, X, X, X, X, X, X, X, X},	// 0010 1011 +
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0010 1100
	{3, 0, 1, 9, 0, 5, 3, 3, 5, 2, X, X, X, X, X, X, X, X, X},	// 0010 1101 +
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0010 1110
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0010 1111

	// パターン3 : 内部/外部輪郭
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0011 0000
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0011 0001
	{1, 4, 1, 5, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0011 0010 *
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0011 0011
	{4, 3, 0, 4, 3, 4, 5, 3, 5, 2, 8, 1, 9, X, X, X, X, X, X},	// 0011 0100 +
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0011 0101
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0011 0110
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0011 0111
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0011 1000
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0011 1001
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0011 1010
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0011 1011
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0011 1100
	{3, 3, 0, 4, 3, 4, 5, 3, 5, 2, X, X, X, X, X, X, X, X, X},	// 0011 1101 *
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0011 1110
	{4, 3, 0, 8, 3, 8, 9, 3, 9, 2, 4, 1, 5, X, X, X, X, X, X},	// 0011 1111 +

	// パターン4 : 内部輪郭(輪郭始点)
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0100 0000
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0100 0001
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0100 0010
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0100 0011
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0100 0100
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0100 0101
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0100 0110
	{3, 0, 10, 3, 0, 1, 6, 1, 2, 6, X, X, X, X, X, X, X, X, X},	// 0100 0111 +
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0100 1000
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0100 1001
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0100 1010
	{3, 0, 6, 3, 0, 1, 6, 1, 2, 10, X, X, X, X, X, X, X, X, X},// 0100 1011 +
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0100 1100
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0100 1101
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0100 1110
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0100 1111

	// パターン5 : 内部/外部輪郭
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0101 0000
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0101 0001
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0101 0010
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0101 0011
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0101 0100
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0101 0101
	{2, 1, 2, 4, 6, 4, 2, X, X, X, X, X, X, X, X, X, X, X, X},	// 0101 0110 *
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0101 0111
	{4, 0, 4, 3, 6, 3, 4, 1, 2, 8, 10, 8, 2, X, X, X, X, X, X},	// 0101 1000 +
	{2, 0, 4, 3, 6, 3, 4, X, X, X, X, X, X, X, X, X, X, X, X},	// 0101 1001 *
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0101 1010
	{4, 0, 8, 3, 10, 3, 8, 1, 2, 4, 6, 4, 2, X, X, X, X, X, X},	// 0101 1011 +
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0101 1100
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0101 1101
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0101 1110
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0101 1111

	// パターン6 : 内部/外部輪郭
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0110 0000
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0110 0001
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0110 0010
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0110 0011
	{1, 2, 6, 5, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0110 0100 *
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0110 0101
	{4, 2, 10, 9, 0, 1, 5, 0, 5, 6, 0, 6, 3, X, X, X, X, X, X},	// 0110 0110 +
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0110 0111
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0110 1000
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0110 1001
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0110 1010
	{3, 0, 1, 5, 0, 5, 6, 0, 6, 3, X, X, X, X, X, X, X, X, X},	// 0110 1011 *
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0110 1100
	{4, 2, 6, 5, 0, 1, 9, 0, 9, 10, 0, 10, 3, X, X, X, X, X, X},// 0110 1101 +
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0110 1110
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 0110 1111

	// パターン7 : 内部/外部輪郭(特殊)
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},		// 0111 0000
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},		// 0111 0001
	{4, 2, 10, 9, 4, 5, 6, 0, 4, 6, 0, 6, 3, X, X, X, X, X, X},		// 0111 0010 *+
	{4, 2, 6, 5, 4, 9, 10, 0, 4, 10, 0, 10, 3, X, X, X, X, X, X},	// 0111 0011 *+
	{4, 1, 5, 4, 8, 9, 6, 0, 8, 6, 0, 6, 3, X, X, X, X, X, X},		// 0111 0100 *+
	{4, 1, 8, 9, 4, 5, 6, 0, 4, 6, 0, 6, 3, X, X, X, X, X, X},		// 0111 0101 *+
	{2, 1, 9, 4, 2, 6, 5, X, X, X, X, X, X, X, X, X, X, X, X},		// 0111 0110 *+
	{2, 1, 5, 4, 2, 6, 9, X, X, X, X, X, X, X, X, X, X, X, X},		// 0111 0111 *+
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},		// 0111 1000
	{5, 4, 5, 6, 0, 4, 6, 0, 6, 3, 1, 9, 8, 2, 10, 12, X, X, X},	// 0111 1001 +
	{5, 2, 6, 5, 4, 9, 10, 0, 4, 10, 0, 10, 3, 1, 12, 8, X, X, X},	// 0111 1010 +
	{5, 4, 5, 6, 0, 4, 6, 0, 6, 3, 2, 10, 9, 1, 12, 8, X, X, X},	// 0111 1011 +
	{5, 1, 5, 4, 2, 6, 9, 8, 12, 10, 0, 8, 10, 0, 10, 3, X, X, X},	// 0111 1100 +
	{5, 1, 5, 4, 8, 9, 6, 0, 8, 6, 0, 6, 3, 2, 10, 12, X, X, X},	// 0111 1101 +
	{5, 2, 6, 5, 1, 9, 4, 8, 12, 10, 0, 8, 10, 0, 10, 3, X, X, X},	// 0111 1110 +
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},		// 0111 1111

	// パターン8 : 内部輪郭(輪郭始点)
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1000 0000
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1000 0001
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1000 0010
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1000 0011
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1000 0100
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1000 0101
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1000 0110
	{3, 0, 1, 7, 1, 2, 7, 2, 3, 11, X, X, X, X, X, X, X, X, X},// 1000 0111 +
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1000 1000
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1000 1001
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1000 1010
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1000 1011
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1000 1100
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1000 1101
	{3, 0, 1, 11, 1, 2, 7, 2, 3, 7, X, X, X, X, X, X, X, X, X},	// 1000 1110 +
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1000 1111

	// パターン9 : 内部/外部輪郭
	{4, 0, 4, 7, 2, 8, 1, 2, 11, 8, 2, 3, 11, X, X, X, X, X, X},// 1001 0000 +
	{1, 0, 4, 7, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1001 0001 *
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1001 0010
	{4, 0, 8, 11, 2, 4, 1, 2, 7, 4, 2, 3, 7, X, X, X, X, X, X},	// 1001 0011 +
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1001 0100
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1001 0101
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1001 0110
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1001 0111
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1001 1000
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1001 1001
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1001 1010
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1001 1011
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1001 1100
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1001 1101
	{3, 2, 4, 1, 2, 7, 4, 2, 3, 7, X, X, X, X, X, X, X, X, X},	// 1001 1110 *
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1001 1111

	// パターン10 : 内部/外部輪郭
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1010 0000
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1010 0001
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1010 0010
	{2, 0, 1, 5, 0, 5, 7, X, X, X, X, X, X, X, X, X, X, X, X},	// 1010 0011 *
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1010 0100
	{4, 0, 1, 9, 0, 9, 11, 3, 7, 5, 3, 5, 2, X, X, X, X, X, X},	// 1010 0101 +
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1010 0110
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1010 0111
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1010 1000
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1010 1001
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1010 1010
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1010 1011
	{2, 3, 7, 5, 3, 5, 2, X, X, X, X, X, X, X, X, X, X, X, X},	// 1010 1100 *
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1010 1101
	{4, 0, 1, 5, 0, 5, 7, 3, 11, 9, 3, 9, 2, X, X, X, X, X, X},	// 1010 1110 +
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},	// 1010 1111

	// パターン11 : 内部/外部輪郭(特殊)
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},		// 1011 0000
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},		// 1011 0001
	{4, 1, 9, 8, 4, 5, 7, 3, 7, 5, 2, 3, 5, X, X, X, X, X, X},		// 1011 0010 *+
	{4, 1, 5, 4, 8, 9, 7, 3, 7, 9, 2, 3, 9, X, X, X, X, X, X},		// 1011 0011 *+
	{4, 0, 4, 7, 8, 5, 11, 3, 11, 5, 2, 3, 5, X, X, X, X, X, X},	// 1011 0100 *+
	{4, 0, 8, 11, 4, 5, 7, 3, 7, 5, 2, 3, 5, X, X, X, X, X, X},		// 1011 0101 *+
	{2, 0, 8, 7, 1, 5, 4, X, X, X, X, X, X, X, X, X, X, X, X},		// 1011 0110 *+
	{2, 0, 4, 7, 1, 5, 8, X, X, X, X, X, X, X, X, X, X, X, X},		// 1011 0111 *+
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},		// 1011 1000
	{5, 4, 5, 7, 3, 7, 5, 2, 3, 5, 0, 8, 11, 1, 9, 12, X, X, X},	// 1011 1001 +
	{5, 1, 5, 4, 8, 9, 7, 3, 7, 9, 2, 3, 9, 0, 12, 11, X, X, X},	// 1011 1010 +
	{5, 4, 5, 7, 3, 7, 5, 2, 3, 5, 1, 9, 8, 0, 12, 11, X, X, X},	// 1011 1011 +
	{5, 0, 4, 7, 1, 5, 8, 12, 9, 11, 3, 11, 9, 2, 3, 9, X, X, X},	// 1011 1100 +
	{5, 0, 4, 7, 8, 5, 11, 3, 11, 5, 2, 3, 5, 1, 9, 12, X, X, X},	// 1011 1101 +
	{5, 1, 5, 4, 0, 8, 7, 12, 9, 11, 3, 11, 9, 2, 3, 9, X, X, X},	// 1011 1110 +
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},		// 1011 1111

	// パターン12 : 内部/外部輪郭
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},		// 1100 0000
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},		// 1100 0001
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},		// 1100 0010
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},		// 1100 0011
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},		// 1100 0100
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},		// 1100 0101
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},		// 1100 0110
	{3, 1, 7, 0, 1, 6, 7, 1, 2, 6, X, X, X, X, X, X, X, X, X},		// 1100 0111 *
	{1, 3, 7, 6, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},		// 1100 1000 *
	{4, 1, 11, 0, 1, 10, 11, 1, 2, 10, 3, 7, 6, X, X, X, X, X, X},	// 1100 1001 +
	{4, 1, 7, 0, 1, 6, 7, 1, 2, 6, 3, 11, 10, X, X, X, X, X, X},	// 1100 1010 +
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},		// 1100 1011
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},		// 1100 1100
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},		// 1100 1101
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},		// 1100 1110
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},		// 1100 1111

	// パターン13 : 内部/外部輪郭(特殊)
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},		// 1101 0000
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},		// 1101 0001
	{4, 0, 8, 11, 4, 6, 7, 2, 6, 4, 1, 2, 4, X, X, X, X, X, X},		// 1101 0010 *+
	{4, 0, 4, 7, 8, 6, 11, 2, 6, 8, 1, 2, 8, X, X, X, X, X, X},		// 1101 0011 *+
	{4, 3, 7, 6, 4, 10, 11, 2, 10, 4, 1, 2, 4, X, X, X, X, X, X},	// 1101 0100 *+
	{4, 3, 11, 10, 4, 6, 7, 2, 6, 4, 1, 2, 4, X, X, X, X, X, X},	// 1101 0101 *+
	{2, 0, 4, 7, 3, 11, 6, X, X, X, X, X, X, X, X, X, X, X, X},		// 1101 0110 *+
	{2, 0, 4, 11, 3, 7, 6, X, X, X, X, X, X, X, X, X, X, X, X},		// 1101 0111 *+
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},		// 1101 1000
	{5, 4, 6, 7, 2, 6, 4, 1, 2, 4, 3, 11, 10, 0, 8, 12, X, X, X},	// 1101 1001 +
	{5, 0, 4, 7, 8, 6, 11, 2, 6, 8, 1, 2, 8, 3, 12, 10, X, X, X},	// 1101 1010 +
	{5, 4, 6, 7, 2, 6, 4, 1, 2, 4, 0, 8, 11, 3, 12, 10, X, X, X},	// 1101 1011 +
	{5, 3, 7, 6, 0, 4, 11, 8, 10, 12, 2, 10, 8, 1, 2, 8, X, X, X},	// 1101 1100 +
	{5, 3, 7, 6, 4, 10, 11, 2, 10, 4, 1, 2, 4, 0, 8, 12, X, X, X},	// 1101 1101 +
	{5, 0, 4, 7, 3, 11, 6, 8, 10, 12, 2, 10, 8, 1, 2, 8, X, X, X},	// 1101 1110 +
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},		// 1101 1111

	// パターン14 : 内部/外部輪郭(特殊)
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},		// 1110 0000
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},		// 1110 0001
	{4, 3, 11, 10, 5, 6, 7, 1, 5, 7, 0, 1, 7, X, X, X, X, X, X},	// 1110 0010 *+
	{4, 3, 7, 6, 5, 10, 11, 1, 5, 11, 0, 1, 11, X, X, X, X, X, X},	// 1110 0011 *+
	{4, 2, 6, 5, 9, 10, 7, 1, 9, 7, 0, 1, 7, X, X, X, X, X, X},		// 1110 0100 *+
	{4, 2, 10, 9, 5, 6, 7, 1, 5, 7, 0, 1, 7, X, X, X, X, X, X},		// 1110 0101 *+
	{2, 3, 7, 6, 2, 10, 5, X, X, X, X, X, X, X, X, X, X, X, X},		// 1110 0110 *+
	{2, 3, 7, 10, 2, 6, 5, X, X, X, X, X, X, X, X, X, X, X, X},		// 1110 0111 *+
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},		// 1110 1000
	{5, 5, 6, 7, 1, 5, 7, 0, 1, 7, 2, 10, 9, 3, 11, 12, X, X, X},	// 1110 1001 +
	{5, 3, 7, 6, 5, 10, 11, 1, 5, 11, 0, 1, 11, 2, 12, 9, X, X, X},	// 1110 1010 +
	{5, 5, 6, 7, 1, 5, 7, 0, 1, 7, 3, 11, 10, 2, 12, 9, X, X, X},	// 1110 1011 +
	{5, 2, 6, 5, 3, 7, 10, 9, 12, 11, 1, 9, 11, 0, 1, 11, X, X, X},	// 1110 1100 +
	{5, 2, 6, 5, 9, 10, 7, 1, 9, 7, 0, 1, 7, 3, 11, 12, X, X, X},	// 1110 1101 +
	{5, 3, 7, 6, 2, 10, 5, 9, 12, 11, 1, 9, 11, 0, 1, 11, X, X, X},	// 1110 1110 +
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},		// 1110 1111

	// パターン15 : 内部/外部輪郭(特殊)
	{2, 0, 4, 7, 2, 6, 5, X, X, X, X, X, X, X, X, X, X, X, X},			// 1111 0000 *
	{3, 3, 7, 6, 1, 5, 4, 2, 10, 9, X, X, X, X, X, X, X, X, X},			// 1111 0001 *
	{3, 2, 6, 5, 3, 7, 10, 1, 9, 4, X, X, X, X, X, X, X, X, X},			// 1111 0010 *
	{3, 3, 7, 6, 2, 10, 5, 1, 9, 4, X, X, X, X, X, X, X, X, X},			// 1111 0011 *
	{3, 1, 5, 4, 2, 6, 9, 3, 7, 10, X, X, X, X, X, X, X, X, X},			// 1111 0100 *
	{3, 1, 5, 4, 3, 7, 6, 2, 10, 9, X, X, X, X, X, X, X, X, X},			// 1111 0101 *
	{3, 2, 6, 5, 1, 9, 4, 3, 7, 10, X, X, X, X, X, X, X, X, X},			// 1111 0110 *
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},			// 1111 0111
	{0, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},			// 1111 1000
	{6, 3, 7, 6, 1, 5, 4, 2, 10, 9, 0, 8, 11, 8, 12, 11, 12, 13, 11},	// 1111 1001 +
	{6, 2, 6, 5, 3, 7, 10, 1, 9, 4, 0, 8, 11, 8, 12, 11, 12, 13, 11},	// 1111 1010 +
	{6, 3, 7, 6, 2, 10, 5, 1, 9, 4, 0, 8, 11, 8, 12, 11, 12, 13, 11},	// 1111 1011 +
	{6, 1, 5, 4, 2, 6, 9, 3, 7, 10, 0, 8, 11, 8, 12, 11, 12, 13, 11},	// 1111 1100 +
	{6, 1, 5, 4, 3, 7, 6, 2, 10, 9, 0, 8, 11, 8, 12, 11, 12, 13, 11},	// 1111 1101 +
	{6, 2, 6, 5, 1, 9, 4, 3, 7, 10, 0, 8, 11, 8, 12, 11, 12, 13, 11},	// 1111 1110 +
	{2, 1, 5, 4, 3, 7, 6, X, X, X, X, X, X, X, X, X, X, X, X},			// 1111 1111 *
};


/*!
 * 各パターンにおけるエッジ頂点リスト
 */
const int g_EdgeTable[16][4] = {
	{X, X, X, X},	// 0
	{4, X, X, X},	// 1
	{5, X, X, X},	// 2
	{4, 5, X, X},	// 3

	{6, X, X, X},	// 4
	{4, 6, X, X},	// 5
	{5, 6, X, X},	// 6
	{4, 5, 6, X},	// 7

	{7, X, X, X},	// 8
	{7, 4, X, X},	// 9
	{5, 7, X, X},	// 10
	{7, 4, 5, X},	// 11

	{6, 7, X, X},	// 12
	{6, 7, 4, X},	// 13
	{5, 6, 7, X},	// 14
	{4, 5, 6, 7},	// 15
};


/*!
 * パターン15におけるノード頂点,エッジ頂点リスト
 */
const int g_NodeTable[4][8] = {
	{3, 1, 2, 3,  4, 5, 6, 7},	// 1110
	{0, 2, 3, 0,  5, 6, 7, 4},	// 1101
	{1, 3, 0, 1,  6, 7, 4, 5},	// 1011
	{2, 0, 1, 2,  7, 4, 5, 6},	// 0111
};

/*!
 * 頂点回転対応テーブル
 */
const int g_VrtRotTable[4][14] = {
	{0, 1, 2, 3,  4, 5, 6, 7,  8, 9, 10, 11,  12, 13}, 
	{1, 2, 3, 0,  5, 6, 7, 4,  9, 10, 11, 8,  12, 13}, 
	{2, 3, 0, 1,  6, 7, 4, 5,  10, 11, 8, 9,  12, 13}, 
	{3, 0, 1, 2,  7, 4, 5, 6,  11, 8, 9, 10,  12, 13}
};

#undef X


#endif // #ifdef _RX_SSM_TABLES_H_
