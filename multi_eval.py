from math import ceil, log
from numpy.fft import fft, ifft
import numpy as np
from poly_arith import poly_mult, poly_div
	
#Going up the product tree in the multipoint evaluation algorithm 
def create_tree(r):
	num = len(r)
	tree = [[]]
	for i in range(num):
		tree[0].append([1,-r[i]]) 
	i = 1
	j = 2
	while j <= num:
		tree.append([])
		for k in range(num // j):
	  		prod = poly_mult(tree[i-1][k*2], tree[i-1][k*2 + 1])    
	  		tree[i].append(prod) 
		j = j*2
		i = i+1
	return tree 
		
#Going down the product tree 
def mult_eval(poly, points):
	num = len(points)
	if num == 1:
		return poly
	j = 1
	while j < num:
		j = 2*j
	eval_points = []
	if num == j:
		for i in range(num):
	  		eval_points.append(points[i])
	else:
		for i in range(num):
	  		eval_points.append(points[i])
		for i in range(j - num):
	  		eval_points.append(0)
	tree = create_tree(eval_points)
	r0 = poly_div(poly, tree[len(tree)-2][0])
	r1 = poly_div(poly, tree[len(tree)-2][1])
	list1 = mult_eval_recurse(r0, eval_points[:j//2], tree, 1, 0)
	list2 = mult_eval_recurse(r1, eval_points[j//2:], tree, 1, 1)
	return (list1 + list2)[:num] 

#A helper recursion function for going down product tree
def mult_eval_recurse(poly, points, tree, depth, index):
	num = len(points)
	if num == 1:
		return list(poly)
	r0 = poly_div(poly, tree[len(tree)-2-depth][2*index])
	r1 = poly_div(poly, tree[len(tree)-2-depth][2*index+1])
	list1 = mult_eval_recurse(r0, points[:num//2], tree, depth + 1, 2*index)
	list2 = mult_eval_recurse(r1, points[num//2:], tree, depth + 1, 2*index+1)
	return (list1 + list2)[:num] 
   
