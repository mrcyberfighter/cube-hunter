#!/usr/bin/python
# -*- coding: utf-8 -*-

#####################################################################################
#                                                                                   #       
# cube hunter: an shoot-them-all game in a cubic meshed playground where the        #
#              pentahedron player have to eliminate all the cubes.                  #
# Copyright (C) 2014 Brüggemann Eddie                                               #
#                                                                                   #
# This file is part of cube hunter.                                                 #
# cube hunter is free software: you can redistribute it and/or modify               #
# it under the terms of the GNU General Public License as published by              #
# the Free Software Foundation, either version 3 of the License, or                 #
# (at your option) any later version.                                               #
#                                                                                   #
# cube hunter is distributed in the hope that it will be useful,                    #  
# but WITHOUT ANY WARRANTY; without even the implied warranty of                    # 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                      #
# GNU General Public License for more details.                                      #
#                                                                                   #
# You should have received a copy of the GNU General Public License                 #
# along with cube hunter. If not, see <http://www.gnu.org/licenses/>                #
#                                                                                   #
#####################################################################################

from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *

from time import sleep
from random import randint,choice

from math import *

from sys import exit

from Tkinter import Tk

import pygame

from pygame.locals import *

class Matrix(object) :
  def __init__(self,vertex=False) :
    if vertex :
      if not isinstance(vertex,Vertex) :
        raise TypeError(Vertex)
      
      self.x=vertex.wx
      self.y=vertex.wy
      self.z=vertex.wz
    else :
      self.x=1.
      self.y=1.
      self.z=1.
      
    self.matrix=range(0,16)
    self.matrix[0]=1      ; self.matrix[4]=0      ; self.matrix[8]=0       ; self.matrix[12]=0 ;
    self.matrix[1]=0      ; self.matrix[5]=1      ; self.matrix[9]=0       ; self.matrix[13]=0 ;
    self.matrix[2]=0      ; self.matrix[6]=0      ; self.matrix[10]=1      ; self.matrix[14]=0 ;
    self.matrix[3]=0      ; self.matrix[7]=0      ; self.matrix[11]=0      ; self.matrix[15]=1 ;
    
    self._move_matrix=range(0,16)                 
    
  
  def translate(self,vector) :
    if not isinstance(vector,tuple) and not isinstance(vector,list) :
      raise TypeError(tuple,list)
    
    self._move_matrix[0]=1 ; self._move_matrix[4]=0 ; self._move_matrix[8]=0  ; self._move_matrix[12]=vector[0]  ;
    self._move_matrix[1]=0 ; self._move_matrix[5]=1 ; self._move_matrix[9]=0  ; self._move_matrix[13]=vector[1]  ;
    self._move_matrix[2]=0 ; self._move_matrix[6]=0 ; self._move_matrix[10]=1 ; self._move_matrix[14]=vector[2]  ;
    self._move_matrix[3]=0 ; self._move_matrix[7]=0 ; self._move_matrix[11]=0 ; self._move_matrix[15]=1          ;
                      
    self._multiply()
    
  def rotate_x(self,degrees) :
    c=cos(radians(degrees))
    s=sin(radians(degrees))
    self._move_matrix[0]=1 ; self._move_matrix[4]=0 ; self._move_matrix[8]=0  ; self._move_matrix[12]=0 ;
    self._move_matrix[1]=0 ; self._move_matrix[5]=c ; self._move_matrix[9]=-s ; self._move_matrix[13]=0 ;
    self._move_matrix[2]=0 ; self._move_matrix[6]=s ; self._move_matrix[10]=c ; self._move_matrix[14]=0 ;
    self._move_matrix[3]=0 ; self._move_matrix[7]=0 ; self._move_matrix[11]=0 ; self._move_matrix[15]=1 ;
    self._multiply()
  
  def rotate_y(self,degrees) :
    c=cos(radians(degrees))
    s=sin(radians(degrees))
    self._move_matrix[0]=c  ; self._move_matrix[4]=0 ; self._move_matrix[8]=s  ; self._move_matrix[12]=0 ;
    self._move_matrix[1]=0  ; self._move_matrix[5]=1 ; self._move_matrix[9]=0  ; self._move_matrix[13]=0 ;
    self._move_matrix[2]=-s ; self._move_matrix[6]=0 ; self._move_matrix[10]=c ; self._move_matrix[14]=0 ;
    self._move_matrix[3]=0  ; self._move_matrix[7]=0 ; self._move_matrix[11]=0 ; self._move_matrix[15]=1 ;
    self._multiply()
  
  
  def rotate_z(self,degrees) :
    c=cos(radians(degrees))
    s=sin(radians(degrees))
    self._move_matrix[0]=c ; self._move_matrix[4]=-s ; self._move_matrix[8]=0  ; self._move_matrix[12]=0 ;
    self._move_matrix[1]=s ; self._move_matrix[5]=c  ; self._move_matrix[9]=0  ; self._move_matrix[13]=0 ;
    self._move_matrix[2]=0 ; self._move_matrix[6]=0  ; self._move_matrix[10]=1 ; self._move_matrix[14]=0 ;
    self._move_matrix[3]=0 ; self._move_matrix[7]=0  ; self._move_matrix[11]=0 ; self._move_matrix[15]=1 ;
    self._multiply()
  
  def load_hardware(self) :
    glMatrixMode(GL_MODELVIEW)
    glLoadMatrixd(self.matrix)
  
  def mult_vertex(self,vertex) :
    if not isinstance(vertex,Vertex) :
      raise TypeError(Vertex)
    
    res_x= vertex.wx * self.matrix[0] + vertex.wy * self.matrix[4] + vertex.wz * self.matrix[8]  + self.matrix[12]
    res_y= vertex.wx * self.matrix[1] + vertex.wy * self.matrix[5] + vertex.wz * self.matrix[9]  + self.matrix[13]
    res_z= vertex.wx * self.matrix[2] + vertex.wy * self.matrix[6] + vertex.wz * self.matrix[10] + self.matrix[14]
    
    return Vertex(res_x,res_y,res_z)   
  
  def mult_vector(self,vector) :
    if not isinstance(vector,Vector) :
      raise TypeError(Vector)
    x= vector.x*self.matrix[0] + vector.y*self.matrix[4] + vector.z*self.matrix[8]
    y= vector.x*self.matrix[1] + vector.y*self.matrix[5] + vector.z*self.matrix[9]
    z= vector.x*self.matrix[2] + vector.y*self.matrix[6] + vector.z*self.matrix[10]
    
    return Vector(x,y,z)
  
  def rotate_vector(self,angle,vector) :
    if not isinstance(vector,Vector) :
      raise TypeError(Vector)
    
    v = vector.mult_vector((1.0 / vector.length),vector)
    x = v.x ; y = v.y ; z = v.z
    
    s = sin(radians(angle))
    c = cos(radians(angle))
    
    self._move_matrix[0]= x*x*(1-c)+c   ; self._move_matrix[4]= x*y*(1-c)-z*s ; self._move_matrix[8] = x*z*(1-c)+y*s ; self._move_matrix[12]=0 ;
    self._move_matrix[1]= x*y*(1-c)+z*s ; self._move_matrix[5]= y*y*(1-c)+c   ; self._move_matrix[9] = y*z*(1-c)-x*s ; self._move_matrix[13]=0 ;
    self._move_matrix[2]= x*z*(1-c)-y*s ; self._move_matrix[6]= y*z*(1-c)+x*s ; self._move_matrix[10]= z*z*(1-c)+c   ; self._move_matrix[14]=0 ;
    self._move_matrix[3]=0              ; self._move_matrix[7]=0              ; self._move_matrix[11]= 0             ; self._move_matrix[15]=1 ;
    self._multiply()
 
  def _multiply(self) :
    x=0
    tmp_matrix=range(0,16)
    while x < 16 :
      value1= x % 4
      value2=(x/4)*4
      tmp_matrix[x]=self._move_matrix[value1] * self.matrix[value2+0] + self._move_matrix[value1+4] * self.matrix[value2+1]+ self._move_matrix[value1+8] * self.matrix[value2+2]+ self._move_matrix[value1+12] *  self.matrix[value2+3]
      x += 1
    
    for v in range(0,16) :
      self.matrix[v]=tmp_matrix[v]

class Vertex(object) :
  def __init__(self,x=False,y=False,z=False,vertexv=False) :
    ''' Object from type Vertex representing an vertice '''
    if vertexv :
      self.wx=vertexv[0]
      self.wy=vertexv[1]
      self.wz=vertexv[2]
    else :
      self.wx=x
      self.wy=y
      self.wz=z
  
  def get_vertex(self) :
    return (self.wx,self.wy,self.wz)     
        
class Vector(object) :
  def __init__(self,x,y,z) :
    ''' Object from type Vector representing an vector '''
    self.x=x
    self.y=y
    self.z=z
    self.length=sqrt(pow(x,2)+pow(y,2)+pow(z,2))
   
  def mult_vector(self,mult,vector) :
    ''' increment vertor length and
        if mult is negativ flip the vector direction. '''
    if not isinstance(vector,Vector) :
      raise TypeError(Vector)
    else :
      return Vector(mult*vector.x,mult*vector.y,mult*vector.z)
  
  
    
class Localview(object) :
  
  def __init__(self,x=0.,y=0.,z=0.) :
    self.pos=Vertex(x,y,z) #instanciation not sure to work as wanted ??? 
    self.right=Vector(1.,0.,0.)
    self.up=Vector(0.,1.,0.)
    self.sight=Vector(0.,0.,1.)
    
  def mult_matrix(self,matrix,local_view) :
    if not isinstance(matrix,Matrix) :
      raise TypeError(Matrix)
    
    local_view.pos=matrix.mult_vertex(local_view.pos)
    local_view.up=matrix.mult_vector(local_view.up)
    local_view.right=matrix.mult_vector(local_view.right)
    local_view.sight=matrix.mult_vector(local_view.sight)
    return local_view
  



def set_point_size(size,is_round=False) :
  if is_round :
    glEnable(GL_POINT_SMOOTH)
    glEnable(GL_BLEND)
    glHint(GL_POINT_SMOOTH_HINT,GL_NICEST)
  else :
    glDisable(GL_POINT_SMOOTH)
    glDisable(GL_BLEND)
    
  glPointSize(size)  



def set_polygon_mode(side=GL_FRONT_AND_BACK,mode=GL_FILL) :
  '''side= GL_FRONT_AND_BACK | GL_FRONT | GL_BACK 
     mode= GL_FILL | GL_LINE | GL_POINT'''
     
  glPolygonMode(side,mode)
  

def set_lined(line_width) :
  glPolygonMode(GL_FRONT_AND_BACK,GL_LINE)
  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA)
  glEnable(GL_LINE_SMOOTH)
  glEnable(GL_BLEND)
  glHint(GL_LINE_SMOOTH_HINT,GL_NICEST)
  glLineWidth(line_width)



def generate_polygon_on_xy(edges,radius,center,offset=0) :
  ''' Return an polygon on surface xy:
           | 
       ____|____
           |
           | 
           
      from edges sides, from radius radius, with offset offset'''
  
  res=[]
  scale=360./edges
  i=0
  while i < edges :
    res.append(Vertex(radius*cos(radians((scale*i)+offset))+center.wx,radius*sin(radians((scale*i)+offset))+center.wy,center.wz))
    
    i += 1
  return res


def gen_player(side_length) :
  ''' Return an array of triangles '''
  
  vertex_list=[Vertex(-side_length/2., -side_length/2., -side_length/2.), 
               Vertex(side_length/2.,  -side_length/2., -side_length/2.), 
               Vertex(side_length/2.,   side_length/2., -side_length/2.), 
               Vertex(-side_length/2.,  side_length/2., -side_length/2.),
               Vertex(-side_length/2., -side_length/2., side_length/2.), 
               Vertex(side_length/2.,  -side_length/2., side_length/2.), 
               Vertex(side_length/2.,   side_length/2., side_length/2.), 
               Vertex(-side_length/2.,  side_length/2., side_length/2.)]

  
  middle=get_middle(vertex_list[0],vertex_list[2])
  distance=get_distance(middle,Vertex(0.,0.,0.))
  
  pointed=Vertex(0.,0.,distance*2.)
  
  
  points=[]
  triangles=[]
  points.append(pointed)
  
  
  
  for v in xrange(0,4) :
    points.append(vertex_list[v])
  
  i=0
  while i < len(vertex_list) :
    tmp=[]
    if i+1 != len(vertex_list) : 
      tmp.append(vertex_list[i])
      tmp.append(vertex_list[i+1])
      tmp.append(pointed)
    else :
      tmp.append(vertex_list[i])
      tmp.append(vertex_list[0])
      tmp.append(pointed)  
    triangles.append(tmp)   
    i += 1 
  
  
  return triangles,points,vertex_list 

def get_distance_as_vertex(vertex1,vertex2) :
  ''' Return the local distance between vertex1 and vertex2. 
      The distance from vertex1 to vertex2 as an Vertex.
  '''    
  if not isinstance(vertex1,Vertex) or not isinstance(vertex2,Vertex) :
    raise TypeError(Vertex)
  
  dist_x=abs(sqrt(pow(vertex1.wx-vertex2.wx,2)))
  dist_y=abs(sqrt(pow(vertex1.wy-vertex2.wy,2)))
  dist_z=abs(sqrt(pow(vertex1.wz-vertex2.wz,2)))
  
  if vertex1.wx > vertex2.wx :
    res_x=-dist_x
  elif vertex1.wx < vertex2.wx :
    res_x=dist_x
  else :
    res_x=0.0
    
  if vertex1.wy > vertex2.wy :
    res_y=-dist_y
  elif vertex1.wy < vertex2.wy :
    res_y=dist_y
  else :
    res_y=0.0
    
  if vertex1.wz > vertex2.wz :
    res_z=-dist_z
  elif vertex1.wz < vertex2.wz :
    res_z=dist_z
  else :
    res_z=0.0
    
  return Vertex(res_x,res_y,res_z)  


def compute_middle_point(points) :
  center=[(max(map(lambda e :e.wx,points))+min(map(lambda e :e.wx,points)))/2.,(max(map(lambda e :e.wy,points))+min(map(lambda e :e.wy,points)))/2.,(max(map(lambda e :e.wz,points))+min(map(lambda e :e.wz,points)))/2.] 
  
  return Vertex(center[0],center[1],center[2])



def get_middle(vertex1,vertex2) :
  if not isinstance(vertex1,Vertex) or not isinstance(vertex1,Vertex) :
    raise TypeError(Vertex)
  
  return Vertex((vertex1.wx+vertex2.wx)/2.,(vertex1.wy+vertex2.wy)/2.,(vertex1.wz+vertex2.wz)/2.)  

def get_distance(vertex1,vertex2) :
  if not isinstance(vertex1,Vertex) or not isinstance(vertex1,Vertex) :
    raise TypeError(Vertex)
  
  return sqrt(pow(vertex1.wx-vertex2.wx,2)+pow(vertex1.wy-vertex2.wy,2)+pow(vertex1.wz-vertex2.wz,2))


class Main_Cube_Base_Class() :
  
  # Base class contains positions inherit from the Main_Cube class.
  
  def initialize_shared_datas(self) :
    self.vertex_list=[Vertex(-32.5/2.,-32.5/2., -32.5/2.), Vertex(32.5/2.,-32.5/2., -32.5/2.), Vertex(32.5/2., 32.5/2., -32.5/2.), Vertex(-32.5/2., 32.5/2., -32.5/2.),
                      Vertex(-32.5/2.,-32.5/2.,  32.5/2.), Vertex(32.5/2.,-32.5/2.,  32.5/2.), Vertex(32.5/2., 32.5/2.,  32.5/2.), Vertex(-32.5/2., 32.5/2.,  32.5/2.)]
    
    
    self.front_left_down  = self.vertex_list[0]
    self.back_left_down   = self.vertex_list[1]
    self.back_right_down  = self.vertex_list[2]
    self.front_right_down = self.vertex_list[3]
    
    self.front_left_up    = self.vertex_list[4]
    self.back_left_up     = self.vertex_list[5]
    self.back_right_up    = self.vertex_list[6]
    self.front_right_up   = self.vertex_list[7]
    
    self.front = [self.vertex_list[0],self.vertex_list[1],self.vertex_list[2],self.vertex_list[3]]
    self.back  = [self.vertex_list[4],self.vertex_list[5],self.vertex_list[6],self.vertex_list[7]]
    
    self.left  = [self.vertex_list[7],self.vertex_list[4],self.vertex_list[0],self.vertex_list[3]]
    self.right = [self.vertex_list[5],self.vertex_list[6],self.vertex_list[2],self.vertex_list[1]]
    
    self.up    = [self.vertex_list[1],self.vertex_list[0],self.vertex_list[4],self.vertex_list[5]]
    self.down  = [self.vertex_list[6],self.vertex_list[2],self.vertex_list[3],self.vertex_list[7]]
    
    # Middles Front-Back-Left between Up and Down.
    self.middle_front_back_at_left_down  = get_middle(self.front_left_down,self.back_left_down)
    self.middle_front_back_at_left_up    = get_middle(self.front_left_up,self.back_left_up)  
    
    self.center_left_side = get_middle(self.middle_front_back_at_left_down,self.middle_front_back_at_left_up)    # Left side middle point.
    
    # Middles Front-Back-Right between Up and Down.
    self.middle_front_back_at_right_down = get_middle(self.front_right_down,self.back_right_down)
    self.middle_front_back_at_right_up   = get_middle(self.front_right_up,self.back_right_up)
    
    self.center_right_side = get_middle(self.middle_front_back_at_right_down,self.middle_front_back_at_right_up) # Right side middle point.
    
    # Middles Front between Left and Right Upside and Downside.
    self.middle_left_right_at_front_down = get_middle(self.front_left_down,self.front_right_down)
    self.middle_left_right_at_front_up   = get_middle(self.front_left_up,self.front_right_up)      
                                                                                                               
    self.center_front_side = get_middle(self.middle_left_right_at_front_down,self.middle_left_right_at_front_up) # Front side middle point.
    
    # Middles Back between Left and Right Upside and Downside.
    self.middle_left_right_at_back_down = get_middle(self.back_left_down,self.back_right_down)
    self.middle_left_right_at_back_up = get_middle(self.back_left_up,self.back_right_up)
    
    self.center_back_side = get_middle(self.middle_left_right_at_back_down,self.middle_left_right_at_back_up)    # Back side middle point.
    
    # Up and Down sides centers.
    self.center_down_side = get_middle(self.middle_front_back_at_left_down,self.middle_front_back_at_right_down) # Down side middle point.
    self.center_up_side   = get_middle(self.middle_front_back_at_left_up,self.middle_front_back_at_right_up)     # Up side middle point.
    
    # Up and Down sides centers.
    self.cube_center = get_middle(self.center_down_side,self.center_up_side)                                     # Cube center, middle point.
    
    # Middles from Up to Down at the 4 corners.
    self.middle_up_down_at_front_right = get_middle(self.front_right_down,self.front_right_up)
    self.middle_up_down_at_front_left  = get_middle(self.front_left_down,self.front_left_up)   
    self.middle_up_down_at_back_right = get_middle(self.back_right_down,self.back_right_up)
    self.middle_up_down_at_back_left  = get_middle(self.back_left_down,self.back_left_up)  
    

class Base_Class() :
  
  # Base class contains positions and moving datas inherit from Cube_IA & Player class.
  
  def initialize_shared_datas(self) :
    self.vertex_list=[Vertex(-32.5/2.,-32.5/2., -32.5/2.), Vertex(32.5/2.,-32.5/2., -32.5/2.), Vertex(32.5/2., 32.5/2., -32.5/2.), Vertex(-32.5/2., 32.5/2., -32.5/2.),
                      Vertex(-32.5/2.,-32.5/2.,  32.5/2.), Vertex(32.5/2.,-32.5/2.,  32.5/2.), Vertex(32.5/2., 32.5/2.,  32.5/2.), Vertex(-32.5/2., 32.5/2.,  32.5/2.)]
    
    
    self.front_left_down  = self.vertex_list[0]
    self.back_left_down   = self.vertex_list[1]
    self.back_right_down  = self.vertex_list[2]
    self.front_right_down = self.vertex_list[3]
    
    self.front_left_up    = self.vertex_list[4]
    self.back_left_up     = self.vertex_list[5]
    self.back_right_up    = self.vertex_list[6]
    self.front_right_up   = self.vertex_list[7]
    
    self.front = [self.vertex_list[0],self.vertex_list[1],self.vertex_list[2],self.vertex_list[3]]
    self.back  = [self.vertex_list[4],self.vertex_list[5],self.vertex_list[6],self.vertex_list[7]]
    
    self.left  = [self.vertex_list[7],self.vertex_list[4],self.vertex_list[0],self.vertex_list[3]]
    self.right = [self.vertex_list[5],self.vertex_list[6],self.vertex_list[2],self.vertex_list[1]]
    
    self.up    = [self.vertex_list[1],self.vertex_list[0],self.vertex_list[4],self.vertex_list[5]]
    self.down  = [self.vertex_list[6],self.vertex_list[2],self.vertex_list[3],self.vertex_list[7]]
    
    # Middles Front-Back-Left between Up and Down.
    self.middle_front_back_at_left_down  = get_middle(self.front_left_down,self.back_left_down)
    self.middle_front_back_at_left_up    = get_middle(self.front_left_up,self.back_left_up)  
    
    self.center_left_side = get_middle(self.middle_front_back_at_left_down,self.middle_front_back_at_left_up)    # Left side middle point.
    
    # Middles Front-Back-Right between Up and Down.
    self.middle_front_back_at_right_down = get_middle(self.front_right_down,self.back_right_down)
    self.middle_front_back_at_right_up   = get_middle(self.front_right_up,self.back_right_up)
    
    self.center_right_side = get_middle(self.middle_front_back_at_right_down,self.middle_front_back_at_right_up) # Right side middle point.
    
    # Middles Front between Left and Right Upside and Downside.
    self.middle_left_right_at_front_down = get_middle(self.front_left_down,self.front_right_down)
    self.middle_left_right_at_front_up   = get_middle(self.front_left_up,self.front_right_up)      
                                                                                                               
    self.center_front_side = get_middle(self.middle_left_right_at_front_down,self.middle_left_right_at_front_up) # Front side middle point.
    
    # Middles Back between Left and Right Upside and Downside.
    self.middle_left_right_at_back_down = get_middle(self.back_left_down,self.back_right_down)
    self.middle_left_right_at_back_up = get_middle(self.back_left_up,self.back_right_up)
    
    self.center_back_side = get_middle(self.middle_left_right_at_back_down,self.middle_left_right_at_back_up)    # Back side middle point.
    
    # Up and Down sides centers.
    self.center_down_side = get_middle(self.middle_front_back_at_left_down,self.middle_front_back_at_right_down) # Down side middle point.
    self.center_up_side   = get_middle(self.middle_front_back_at_left_up,self.middle_front_back_at_right_up)     # Up side middle point.
    
    # Up and Down sides centers.
    self.cube_center = get_middle(self.center_down_side,self.center_up_side)                                     # Cube center, middle point.
    
    # Middles from Up to Down at the 4 corners.
    self.middle_up_down_at_front_right = get_middle(self.front_right_down,self.front_right_up)
    self.middle_up_down_at_front_left  = get_middle(self.front_left_down,self.front_left_up)   
    self.middle_up_down_at_back_right = get_middle(self.back_right_down,self.back_right_up)
    self.middle_up_down_at_back_left  = get_middle(self.back_left_down,self.back_left_up)  
    
    self.load_move_computing_datas()
    self.load_move_positions_datas()
    
  def load_move_computing_datas(self) :  
    
    ''' Every attribe is an dict contains:
        The self position                 key : 'pos'                         as value: the position in the main cube and an id (index in the dicts conataining list).
        The moveto orientation position   key : 'F' (Forward), 'B' (Back)...  as value: the slide line end position in the main cube and an id (index in the dicts conataining list).
    '''    
    #####################################################################################
    # Sides centers points.
    self.pos_cube_center = { 'pos' : (self.cube_center,13),                               
                             'F' : (self.center_back_side,22)  , 
                             'B' : (self.center_front_side,4),
                             'R' : (self.center_right_side,14) , 
                             'L' : (self.center_left_side,12),
                             'U' : (self.center_up_side,10)    , 
                             'D' : (self.center_down_side,16) } 
    
    self.pos_back_center = { 'pos' : (self.center_back_side,22),
                             'F' : None, 
                             'B' : (self.cube_center,13),
                             'R' : (self.middle_up_down_at_back_right,23),
                             'L' : (self.middle_up_down_at_back_left,21),
                             'U' : (self.middle_left_right_at_back_up,19),
                             'D' : (self.middle_left_right_at_back_down,25) }
    
    self.pos_front_center = { 'pos' : (self.center_front_side,4),
                              'F' : (self.cube_center,13),
                              'B' : None,
                              'R' : (self.middle_up_down_at_front_right,5),
                              'L' : (self.middle_up_down_at_front_left,3),
                              'U' : (self.middle_left_right_at_front_up,1),
                              'D' : (self.middle_left_right_at_front_down,7) }
    
    self.pos_left_center = { 'pos' : (self.center_left_side,12),
                              'F'  : (self.middle_up_down_at_back_left,21),
                              'B'  : (self.middle_up_down_at_front_left,3),
                              'R'  : (self.cube_center,13),
                              'L'  : None,
                              'U'  : (self.middle_front_back_at_left_up,9),
                              'D'  : (self.middle_front_back_at_left_down,15) }
    self.pos_right_center = { 'pos' : (self.center_right_side,14),              
                              'F' : (self.middle_up_down_at_back_right,23),
                              'B' : (self.middle_up_down_at_front_right,5),
                              'R' : None,
                              'L' : (self.cube_center,13),
                              'U' : (self.middle_front_back_at_right_up,11),
                              'D' : (self.middle_front_back_at_right_down,17)  }
    
    self.pos_up_center = { 'pos' : (self.center_up_side,10),
                           'F' : (self.middle_left_right_at_back_up,19),
                           'B' : (self.middle_left_right_at_front_up,1),
                           'R' : (self.middle_front_back_at_right_up,11),
                           'L' : (self.middle_front_back_at_left_up,9),
                           'U' : None,
                           'D' : (self.cube_center,13)                    }
    
    self.pos_down_center = { 'pos' : (self.center_down_side,16),
                             'F' : (self.middle_left_right_at_back_down,25),
                             'B' : (self.middle_left_right_at_front_down,7),
                             'R' : (self.middle_front_back_at_right_down,17),
                             'L' : (self.middle_front_back_at_left_down,15),
                             'U' : (self.cube_center,13),
                             'D':  None                               }
    
    ############################################################################
    ############################################################################
    # Down middles points on front,back,left and right side.
    self.pos_front_down_middle = { 'pos' : (self.middle_left_right_at_front_down,7),
                                   'F' : (self.center_down_side,16),
                                   'B' : None,
                                   'R' : (self.front_right_down,8),
                                   'L' : (self.front_left_down,6),
                                   'U' : (self.center_front_side,4),
                                   'D': None                          }
    
    self.pos_back_down_middle = { 'pos' : (self.middle_left_right_at_back_down,25),
                                  'F' : None,
                                  'B' : (self.center_down_side,16),
                                  'R' : (self.back_right_down,26),
                                  'L' : (self.back_left_down,24),
                                  'U' : (self.center_back_side,22),
                                  'D' : None                          }
    
    self.pos_left_down_middle = { 'pos' : (self.middle_front_back_at_left_down,15),
                                  'F' : (self.back_left_down,24),
                                  'B' : (self.front_left_down,6),
                                  'R' : (self.center_down_side,16),
                                  'L' : None,
                                  'U' : (self.center_left_side,12),
                                  'D' : None,                        }
    
    self.pos_right_down_middle = { 'pos' : (self.middle_front_back_at_right_down,17),
                                   'F' : (self.back_right_down,26),
                                   'B' : (self.front_right_down,8),
                                   'R' : None,
                                   'L' : (self.center_down_side,16),
                                   'U' : (self.center_right_side,14),
                                   'D' : None                        }
    
    ############################################################################
    ############################################################################
    # Up middles points on front,back,left and right side.
    self.pos_front_up_middle = { 'pos' : (self.middle_left_right_at_front_up,1),
                                 'F' : (self.center_up_side,10),
                                 'B' : None,
                                 'R' : (self.front_right_up,2),
                                 'L' : (self.front_left_up,0),
                                 'U' : None,
                                 'D' : (self.center_front_side,4)       }
    
    self.pos_back_up_middle = { 'pos' : (self.middle_left_right_at_back_up,19),
                                'F' : None,
                                'B' : (self.center_up_side,10),
                                'R' : (self.back_right_up,20),
                                'L' : (self.back_left_up,18),
                                'U' : None,
                                'D' : (self.center_back_side,22)         }
    
    self.pos_left_up_middle = { 'pos' : (self.middle_front_back_at_left_up,9),
                                'F' : (self.back_left_up,18),
                                'B' : (self.front_left_up,0),
                                'R' : (self.center_up_side,10),
                                'L' : None,
                                'U' : None,
                                'D' : (self.center_left_side,12)        }
    
    self.pos_right_up_middle = { 'pos' : (self.middle_front_back_at_right_up,11),
                                 'F' : (self.back_right_up,20),
                                 'B' : (self.front_right_up,2),
                                 'R' : None,
                                 'L' : (self.center_up_side,10),
                                 'U' : None,
                                 'D' : (self.center_right_side,14)      }
    
    ############################################################################
    ############################################################################
    # Front middles points on right and left side.
    self.pos_front_right_middle = { 'pos' : (self.middle_up_down_at_front_right,5),
                                   'F' : (self.center_right_side,14),
                                   'B' : None,
                                   'R' : None,
                                   'L' : (self.center_front_side,4),
                                   'U' : (self.front_right_up,2),
                                   'D' : (self.front_right_down,8)     }
   
    self.pos_front_left_middle = { 'pos' : (self.middle_up_down_at_front_left,3),
                                  'F' : (self.center_left_side,12),
                                  'B' : None,
                                  'R' : (self.center_front_side,4),
                                  'L' : None,
                                  'U' : (self.front_left_up,0),
                                  'D' : (self.front_left_down,6)       }
    
    ############################################################################
    ###########################################################################
    # Back middles points on right and left side.
    self.pos_back_right_middle = { 'pos' : (self.middle_up_down_at_back_right,23),
                                  'F' : None,
                                  'B' : (self.center_right_side,14),
                                  'R' : None,
                                  'L' : (self.center_back_side,22),
                                  'U' : (self.back_right_up,20),
                                  'D' : (self.back_right_down,26)       }
   
    self.pos_back_left_middle = { 'pos' : (self.middle_up_down_at_back_left,21),
                                 'F' : None,
                                 'B' : (self.center_left_side,12),
                                 'R' : (self.center_back_side,22),
                                 'L' : None,
                                 'U' : (self.back_left_up,18),
                                 'D' : (self.back_left_down,24)        }
   
    ###########################################################################
    ###########################################################################
    # Front side cube corners.
   
    self.pos_front_left_up = { 'pos' : (self.front_left_up,0),
                              'F' : (self.middle_front_back_at_left_up,9),
                              'B' : None,
                              'R' : (self.middle_left_right_at_front_up,1),
                              'L' : None,
                              'U' : None,
                              'D' : (self.middle_up_down_at_front_left,3) }
   
    self.pos_front_left_down = { 'pos' : (self.front_left_down,6),         
                                'F' : (self.middle_front_back_at_left_down,15),
                                'B' : None,
                                'R' : (self.middle_left_right_at_front_down,7),
                                'L' : None,
                                'U' : (self.middle_up_down_at_front_left,3),
                                'D': None                             }
   
    self.pos_front_right_up = { 'pos' : (self.front_right_up,2),
                               'F' : (self.middle_front_back_at_right_up,11),
                               'D' : (self.middle_up_down_at_front_right,5),
                               'R' : None,
                               'L' : (self.middle_left_right_at_front_up,1),
                               'U' : None,
                               'D' : (self.middle_up_down_at_front_right,5) }
   
    self.pos_front_right_down = { 'pos' : (self.front_right_down,8),
                                 'F' : (self.middle_front_back_at_right_down,17),
                                 'B' : None,
                                 'R' : None,
                                 'L' : (self.middle_left_right_at_front_down,7),
                                 'U' : (self.middle_up_down_at_front_right,5),
                                 'D' : None                             }
   
    ###########################################################################
    # Back side cube corners.
   
    self.pos_back_left_up = { 'pos' : (self.back_left_up,18),
                             'F' : None,
                             'B' : (self.middle_front_back_at_left_up,9),
                             'R' : (self.middle_left_right_at_back_up,19),
                             'L' : None,
                             'U' : None,
                             'D' : (self.middle_up_down_at_back_left,21)   }
   
    self.pos_back_left_down = { 'pos' : (self.back_left_down,24),
                               'F' : None,
                               'B' : (self.middle_front_back_at_left_down,15),
                               'R' : (self.middle_left_right_at_back_down,25),
                               'L' : None,
                               'U' : (self.middle_up_down_at_back_left,21),
                               'D' : None                              }
   
    self.pos_back_right_up = { 'pos' : (self.back_right_up,20),
                              'F' : None,
                              'B' : (self.middle_front_back_at_right_up,11),
                              'R' : None,
                              'L' : (self.middle_left_right_at_back_up,19),
                              'U' : None,
                              'D' : (self.middle_up_down_at_back_right,23) }
   
    self.pos_back_right_down = { 'pos' : (self.back_right_down,26),
                                'F' : None,
                                'B' : (self.middle_front_back_at_right_down,17),
                                'R' : None,
                                'L' : (self.middle_left_right_at_back_down,25),
                                'U' : (self.middle_up_down_at_back_right,23),
                                'D': None                             }
  
  def load_move_positions_datas(self) : 
    ''' Generate an array for moving contains the dicts define above. '''
    self.position=range(0,27)
   
    self.position[0]=self.pos_front_left_up     ; self.position[1]=self.pos_front_up_middle   ; self.position[2]=self.pos_front_right_up     ;
    self.position[3]=self.pos_front_left_middle ; self.position[4]=self.pos_front_center      ; self.position[5]=self.pos_front_right_middle ;
    self.position[6]=self.pos_front_left_down   ; self.position[7]=self.pos_front_down_middle ; self.position[8]=self.pos_front_right_down   ;
   
    self.position[9]=self.pos_left_up_middle    ; self.position[10]=self.pos_up_center        ; self.position[11]=self.pos_right_up_middle   ;
    self.position[12]=self.pos_left_center      ; self.position[13]=self.pos_cube_center      ; self.position[14]=self.pos_right_center      ;
    self.position[15]=self.pos_left_down_middle ; self.position[16]=self.pos_down_center      ; self.position[17]=self.pos_right_down_middle ;
   
    self.position[18]=self.pos_back_left_up     ; self.position[19]=self.pos_back_up_middle   ; self.position[20]=self.pos_back_right_up     ;
    self.position[21]=self.pos_back_left_middle ; self.position[22]=self.pos_back_center      ; self.position[23]=self.pos_back_right_middle ;
    self.position[24]=self.pos_back_left_down   ; self.position[25]=self.pos_back_down_middle ; self.position[26]=self.pos_back_right_down   ;
    
    # An dict for random choosing end slide positions:
    # key: the current position.
    # values: the index in the array from possible end positions.
    self.start_direction_destination_random={ 0  : [1,2,3,6,9,18] , 
					      1  : [0,2,4,7,10,19],
					      2  : [0,1,5,8,11,20],
					      3  : [0,3,4,5,12,21],
					      4  : [1,3,5,7,13,22],
					      5  : [2,3,4,8,14,23],
					      6  : [0,3,7,8,15,24],
					      7  : [1,4,6,8,16,25],
					      8  : [2,5,6,7,17,26],
					      9  : [0,10,11,12,15,18],
					      10 : [1,9,11,13,16,19],
					      11 : [2,9,10,14,17,20],
					      12 : [3,9,13,14,15,21],
					      13 : [4,10,12,14,16,22],
					      14 : [5,11,12,13,17,23],
					      15 : [6,9,12,16,17,24],
					      16 : [7,10,13,15,17,25],
					      17 : [8,11,14,15,16,26],
					      18 : [0,9,19,20,21,24],
					      19 : [1,10,18,20,22,25],
					      20 : [2,11,18,19,23,26],
					      21 : [3,12,18,22,23,24],
					      22 : [4,13,19,21,23,25],
					      23 : [5,14,20,21,22,26],
					      24 : [6,15,18,21,25,26],
					      25 : [7,16,19,22,24,26],
					      26 : [8,17,20,23,24,25] 
					      }
  
  
  
  def _compute_line_unit(self,vertex1,vertex2) :
    ''' Function to compute every position from the slide lines. '''
    res=[]
    res.append(vertex1)
    res.append(get_middle(vertex1,vertex2))
    res.append(vertex2)
    
    
    limit= int(get_distance(vertex1,vertex2))  
    breaker=0
    while breaker < limit :
      breaker,res=self._computing_line_unit(res)
    
    return res
    
  
  def _computing_line_unit(self,res) :
    ''' Function to compute one position from the slide line and the break of the computing. '''
    tmp_res=res
    res=res[0::]
    i=0
    ii=1
    while i < len(res) :
      if not i == len(res)-1 :
        tmp_res.insert(i+ii,get_middle(res[i],res[i+1]))
        
      else :
        break
      i += 1
      ii += 1
      
    res=tmp_res
    
    return len(tmp_res),res   

class Player_explosing() :
  # Player exploding class inherit from Player class.
  def config_player_explose(self) :
    # Prepare the sounds.
    self.explosion_sound_object=pygame.mixer.Sound("/usr/share/cube-hunter/Sound/Effects/Explosion/Explosion.wav")
    self.hit_sound_object=pygame.mixer.Sound("/usr/share/cube-hunter/Sound/Effects/shoot_hit_cube/Shoot_hit_cube.wav")
    
  def update_to_explose(self) :
    
    # The player 4 quads corners vertex from the pentahedron.
    self.vertex_point_1 = self.player_display[0][0]
    self.vertex_point_2 = self.player_display[1][0]
    self.vertex_point_3 = self.player_display[2][0]
    self.vertex_point_4 = self.player_display[3][0]
    
    
    
    if self.player_direction == "F" :
      
      # Opposite points from the player 4 quads corners vertex from the pentahedron. 
      self.vertex_point_5=Vertex(self.vertex_point_1.wx+1.75,self.vertex_point_1.wy,self.vertex_point_1.wz)
      self.vertex_point_6=Vertex(self.vertex_point_2.wx+1.75,self.vertex_point_2.wy,self.vertex_point_2.wz)
      self.vertex_point_7=Vertex(self.vertex_point_3.wx+1.75,self.vertex_point_3.wy,self.vertex_point_3.wz)
      self.vertex_point_8=Vertex(self.vertex_point_4.wx+1.75,self.vertex_point_4.wy,self.vertex_point_4.wz)
      
      # Prepare explosion lines from the player.
      self.explosion_to_point_1=Vertex(self.vertex_point_1.wx-2,self.vertex_point_1.wy-2,self.vertex_point_1.wz+2)
      self.explosion_to_point_2=Vertex(self.vertex_point_2.wx-2,self.vertex_point_2.wy-2,self.vertex_point_2.wz-2)
      self.explosion_to_point_3=Vertex(self.vertex_point_3.wx-2,self.vertex_point_3.wy+2,self.vertex_point_3.wz-2)
      self.explosion_to_point_4=Vertex(self.vertex_point_4.wx-2,self.vertex_point_4.wy+2,self.vertex_point_4.wz+2)
      
      self.explosion_to_point_5=Vertex(self.vertex_point_1.wx+2,self.vertex_point_1.wy-2,self.vertex_point_1.wz+2)
      self.explosion_to_point_6=Vertex(self.vertex_point_2.wx+2,self.vertex_point_2.wy-2,self.vertex_point_2.wz-2)
      self.explosion_to_point_7=Vertex(self.vertex_point_3.wx+2,self.vertex_point_3.wy+2,self.vertex_point_3.wz-2)
      self.explosion_to_point_8=Vertex(self.vertex_point_4.wx+2,self.vertex_point_4.wy+2,self.vertex_point_4.wz+2)
      
     
    elif self.player_direction == "B" :
      
      # Opposite points from the player 4 quads corners vertex from the pentahedron. 
      self.vertex_point_5=Vertex(self.vertex_point_1.wx-1.75,self.vertex_point_1.wy,self.vertex_point_1.wz)
      self.vertex_point_6=Vertex(self.vertex_point_2.wx-1.75,self.vertex_point_2.wy,self.vertex_point_2.wz)
      self.vertex_point_7=Vertex(self.vertex_point_3.wx-1.75,self.vertex_point_3.wy,self.vertex_point_3.wz)
      self.vertex_point_8=Vertex(self.vertex_point_4.wx-1.75,self.vertex_point_4.wy,self.vertex_point_4.wz)
      
      # Prepare explosion lines from the player.
      self.explosion_to_point_1=Vertex(self.vertex_point_1.wx+2,self.vertex_point_1.wy+2,self.vertex_point_1.wz+2)
      self.explosion_to_point_2=Vertex(self.vertex_point_2.wx+2,self.vertex_point_2.wy+2,self.vertex_point_2.wz-2)
      self.explosion_to_point_3=Vertex(self.vertex_point_3.wx+2,self.vertex_point_3.wy-2,self.vertex_point_3.wz-2)
      self.explosion_to_point_4=Vertex(self.vertex_point_4.wx+2,self.vertex_point_4.wy-2,self.vertex_point_4.wz+2)
      
      self.explosion_to_point_5=Vertex(self.vertex_point_1.wx-2,self.vertex_point_1.wy+2,self.vertex_point_1.wz+2)
      self.explosion_to_point_6=Vertex(self.vertex_point_2.wx-2,self.vertex_point_2.wy+2,self.vertex_point_2.wz-2)
      self.explosion_to_point_7=Vertex(self.vertex_point_3.wx-2,self.vertex_point_3.wy-2,self.vertex_point_3.wz-2)
      self.explosion_to_point_8=Vertex(self.vertex_point_4.wx-2,self.vertex_point_4.wy-2,self.vertex_point_4.wz+2)
      
    elif self.player_direction == "R" :
      
      # Opposite points from the player 4 quads corners vertex from the pentahedron. 
      self.vertex_point_5=Vertex(self.vertex_point_1.wx,self.vertex_point_1.wy+1.75,self.vertex_point_1.wz)
      self.vertex_point_6=Vertex(self.vertex_point_2.wx,self.vertex_point_2.wy+1.75,self.vertex_point_2.wz)
      self.vertex_point_7=Vertex(self.vertex_point_3.wx,self.vertex_point_3.wy+1.75,self.vertex_point_3.wz)
      self.vertex_point_8=Vertex(self.vertex_point_4.wx,self.vertex_point_4.wy+1.75,self.vertex_point_4.wz)
      
      # Prepare explosion lines from the player.
      self.explosion_to_point_1=Vertex(self.vertex_point_1.wx+2,self.vertex_point_1.wy-2,self.vertex_point_1.wz+2)
      self.explosion_to_point_2=Vertex(self.vertex_point_2.wx+2,self.vertex_point_2.wy-2,self.vertex_point_2.wz-2)
      self.explosion_to_point_3=Vertex(self.vertex_point_3.wx-2,self.vertex_point_3.wy-2,self.vertex_point_3.wz-2)
      self.explosion_to_point_4=Vertex(self.vertex_point_4.wx-2,self.vertex_point_4.wy-2,self.vertex_point_4.wz+2)
      
      self.explosion_to_point_5=Vertex(self.vertex_point_1.wx+2,self.vertex_point_1.wy+2,self.vertex_point_1.wz+2)
      self.explosion_to_point_6=Vertex(self.vertex_point_2.wx+2,self.vertex_point_2.wy+2,self.vertex_point_2.wz-2)
      self.explosion_to_point_7=Vertex(self.vertex_point_3.wx-2,self.vertex_point_3.wy+2,self.vertex_point_3.wz-2)
      self.explosion_to_point_8=Vertex(self.vertex_point_4.wx-2,self.vertex_point_4.wy+2,self.vertex_point_4.wz+2) 
     
    elif self.player_direction == "L" :
      
      # Opposite points from the player 4 quads corners vertex from the pentahedron. 
      self.vertex_point_5=Vertex(self.vertex_point_1.wx,self.vertex_point_1.wy-1.75,self.vertex_point_1.wz)
      self.vertex_point_6=Vertex(self.vertex_point_2.wx,self.vertex_point_2.wy-1.75,self.vertex_point_2.wz)
      self.vertex_point_7=Vertex(self.vertex_point_3.wx,self.vertex_point_3.wy-1.75,self.vertex_point_3.wz)
      self.vertex_point_8=Vertex(self.vertex_point_4.wx,self.vertex_point_4.wy-1.75,self.vertex_point_4.wz)
      
      # Prepare explosion lines from the player.
      self.explosion_to_point_1=Vertex(self.vertex_point_1.wx-2,self.vertex_point_1.wy+2,self.vertex_point_1.wz+2)
      self.explosion_to_point_2=Vertex(self.vertex_point_2.wx-2,self.vertex_point_2.wy+2,self.vertex_point_2.wz-2)
      self.explosion_to_point_3=Vertex(self.vertex_point_3.wx+2,self.vertex_point_3.wy+2,self.vertex_point_3.wz-2)
      self.explosion_to_point_4=Vertex(self.vertex_point_4.wx+2,self.vertex_point_4.wy+2,self.vertex_point_4.wz+2)
      
      self.explosion_to_point_5=Vertex(self.vertex_point_1.wx-2,self.vertex_point_1.wy-2,self.vertex_point_1.wz+2)
      self.explosion_to_point_6=Vertex(self.vertex_point_2.wx-2,self.vertex_point_2.wy-2,self.vertex_point_2.wz-2)
      self.explosion_to_point_7=Vertex(self.vertex_point_3.wx+2,self.vertex_point_3.wy-2,self.vertex_point_3.wz-2)
      self.explosion_to_point_8=Vertex(self.vertex_point_4.wx+2,self.vertex_point_4.wy-2,self.vertex_point_4.wz+2) 
      
    elif self.player_direction == "U" :
      
      # Opposite points from the player 4 quads corners vertex from the pentahedron. 
      self.vertex_point_5=Vertex(self.vertex_point_1.wx,self.vertex_point_1.wy,self.vertex_point_1.wz+1.75)
      self.vertex_point_6=Vertex(self.vertex_point_2.wx,self.vertex_point_2.wy,self.vertex_point_2.wz+1.75)
      self.vertex_point_7=Vertex(self.vertex_point_3.wx,self.vertex_point_3.wy,self.vertex_point_3.wz+1.75)
      self.vertex_point_8=Vertex(self.vertex_point_4.wx,self.vertex_point_4.wy,self.vertex_point_4.wz+1.75)
      
      # Prepare explosion lines from the player.
      self.explosion_to_point_1=Vertex(self.vertex_point_1.wx-2,self.vertex_point_1.wy-2,self.vertex_point_1.wz-2)
      self.explosion_to_point_2=Vertex(self.vertex_point_2.wx+2,self.vertex_point_2.wy-2,self.vertex_point_2.wz-2)
      self.explosion_to_point_3=Vertex(self.vertex_point_3.wx+2,self.vertex_point_3.wy+2,self.vertex_point_3.wz-2)
      self.explosion_to_point_4=Vertex(self.vertex_point_4.wx-2,self.vertex_point_4.wy+2,self.vertex_point_4.wz-2)
      
      self.explosion_to_point_5=Vertex(self.vertex_point_1.wx-2,self.vertex_point_1.wy-2,self.vertex_point_1.wz+2)
      self.explosion_to_point_6=Vertex(self.vertex_point_2.wx+2,self.vertex_point_2.wy-2,self.vertex_point_2.wz+2)
      self.explosion_to_point_7=Vertex(self.vertex_point_3.wx+2,self.vertex_point_3.wy+2,self.vertex_point_3.wz+2)
      self.explosion_to_point_8=Vertex(self.vertex_point_4.wx-2,self.vertex_point_4.wy+2,self.vertex_point_4.wz+2) 
      
    elif self.player_direction == "D" :
      
      # Opposite points from the player 4 quads corners vertex from the pentahedron. 
      self.vertex_point_5=Vertex(self.vertex_point_1.wx,self.vertex_point_1.wy,self.vertex_point_1.wz-1.75)
      self.vertex_point_6=Vertex(self.vertex_point_2.wx,self.vertex_point_2.wy,self.vertex_point_2.wz-1.75)
      self.vertex_point_7=Vertex(self.vertex_point_3.wx,self.vertex_point_3.wy,self.vertex_point_3.wz-1.75)
      self.vertex_point_8=Vertex(self.vertex_point_4.wx,self.vertex_point_4.wy,self.vertex_point_4.wz-1.75)
      
      # Prepare explosion lines from the player.
      self.explosion_to_point_1=Vertex(self.vertex_point_1.wx+2,self.vertex_point_1.wy-2,self.vertex_point_1.wz+2)
      self.explosion_to_point_2=Vertex(self.vertex_point_2.wx-2,self.vertex_point_2.wy-2,self.vertex_point_2.wz+2)
      self.explosion_to_point_3=Vertex(self.vertex_point_3.wx-2,self.vertex_point_3.wy+2,self.vertex_point_3.wz+2)
      self.explosion_to_point_4=Vertex(self.vertex_point_4.wx+2,self.vertex_point_4.wy+2,self.vertex_point_4.wz+2)
      
      self.explosion_to_point_5=Vertex(self.vertex_point_1.wx+2,self.vertex_point_1.wy-2,self.vertex_point_1.wz-2)
      self.explosion_to_point_6=Vertex(self.vertex_point_2.wx-2,self.vertex_point_2.wy-2,self.vertex_point_2.wz-2)
      self.explosion_to_point_7=Vertex(self.vertex_point_3.wx-2,self.vertex_point_3.wy+2,self.vertex_point_3.wz-2)
      self.explosion_to_point_8=Vertex(self.vertex_point_4.wx+2,self.vertex_point_4.wy+2,self.vertex_point_4.wz-2)    
     
    
    # Compute explosion lines from the player.
    self.explose_line_1 = self._compute_line_unit(self.vertex_point_1,self.explosion_to_point_1)
    self.explose_line_2 = self._compute_line_unit(self.vertex_point_2,self.explosion_to_point_2)
    self.explose_line_3 = self._compute_line_unit(self.vertex_point_3,self.explosion_to_point_3)
    self.explose_line_4 = self._compute_line_unit(self.vertex_point_4,self.explosion_to_point_4)
    self.explose_line_5 = self._compute_line_unit(self.vertex_point_5,self.explosion_to_point_5)
    self.explose_line_6 = self._compute_line_unit(self.vertex_point_6,self.explosion_to_point_6)
    self.explose_line_7 = self._compute_line_unit(self.vertex_point_7,self.explosion_to_point_7)
    self.explose_line_8 = self._compute_line_unit(self.vertex_point_8,self.explosion_to_point_8)
    
    self.explode_lines=[self.explose_line_1,self.explose_line_2,self.explose_line_3,self.explose_line_4,self.explose_line_5,self.explose_line_6,self.explose_line_7,self.explose_line_8]
  
  def explose(self) :
    ''' Player exploding display function. '''
    
    self.explosion_sound_object.play()
    
    max_len=len(self.explode_lines[0])
    idx=0
    size=self.size*3
    
    while idx < max_len :
      set_point_size(size,False)
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
      # Display the playground.
      main_cube.display()
	  
      for v in self.explode_lines :
	# Display one item per corner explosion line.
        glColor3ubv(self.self_color)
        glBegin(GL_POINTS)
        glVertex3fv(v[idx].get_vertex())
        glEnd()
        sleep(0.02)
	  
      sleep(0.15)
      
      idx += 1
      size -= 1
      
      pygame.display.flip()
  
  def is_ready_to_explose(self) :
    ''' Function to check if the player should explode. '''
    
    self.life_points -= 1
    
    self.hit_sound_object.play()
    
    if self.life_points == 0 :
      
      self.update_to_explose()
      self.explose()
      
      self.lifes -= 1
      
      self.is_dead=True
      
      self.life_points=5
      
      self.lose_a_life()
      
      
      
      if self.lifes == 0 :
	player_lose_text=Text(size=1.05,line_width=6,color=(255,0,0))
        player_lose_text.construct_string("you lose",Vertex(-15.0,-5.0,0.0))
	for v in xrange(0,5) :
	  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT) 
	  main_cube.display()
	  main_cube.manage_inside_room()
	  glPushMatrix()
      
	  camera_view=Localview()
	  m=Matrix()                                  # Build the camera position 
	  m.translate((0.0,0.0,0.0))                  # Build the camera position 
	  m.rotate_vector(0,camera_view.right)        # Build the camera position 
	  m.rotate_vector(0,camera_view.up)           # Build the camera position 
	  m.rotate_vector(0,camera_view.sight)        # Build the camera position 
	  m.load_hardware() 
	  
	  player_lose_text.display()
	  
	  glPopMatrix()
	  
	  pygame.display.flip()    
          sleep(0.15*5)
          
	event_to_send=pygame.event.Event(USEREVENT,{"code":"end_game"})
	pygame.event.post(event_to_send)
	
	
    self.is_hit=False
    sleep(0.15)

class Cube_explosing() :
  # Class for an cube explosion inherit from every cube.
  def config_cube_explose(self) :
    # Register cube explsion sound.
    self.explosion_sound_object=pygame.mixer.Sound("/usr/share/cube-hunter/Sound/Effects/Explosion/Explosion.wav")
    
  def update_to_explose(self) :
    ''' Updating explosion datas in relationship to the cube position. '''
    
    self.front_left_down  = self.self_vertex_list[0]
    self.back_left_down   = self.self_vertex_list[1]
    self.back_right_down  = self.self_vertex_list[2]
    self.front_right_down = self.self_vertex_list[3]
    
    self.front_left_up    = self.self_vertex_list[4]
    self.back_left_up     = self.self_vertex_list[5]
    self.back_right_up    = self.self_vertex_list[6]
    self.front_right_up   = self.self_vertex_list[7]
    
    self.front_left_down_explose_vertex=Vertex(self.front_left_down.wx - 2,self.front_left_down.wy + 2,self.front_left_down.wz - 1)
    self.front_left_down_explose_line = self._compute_line_unit(self.front_left_down,self.front_left_down_explose_vertex)
    
    self.back_left_down_explose_vertex=Vertex(self.back_left_down.wx - 2,self.back_left_down.wy + 2,self.back_left_down.wz + 1)
    self.back_left_down_explose_line = self._compute_line_unit(self.back_left_down,self.back_left_down_explose_vertex)
    
    self.back_right_down_explose_vertex=Vertex(self.back_right_down.wx + 2,self.back_right_down.wy + 2,self.back_right_down.wz + 1)
    self.back_right_down_explose_line = self._compute_line_unit(self.back_right_down,self.back_right_down_explose_vertex)
    
    self.front_right_down_explose_vertex=Vertex(self.front_right_down.wx + 2,self.front_right_down.wy + 2,self.front_right_down.wz - 1)
    self.front_right_down_explose_line = self._compute_line_unit(self.front_right_down,self.front_right_down_explose_vertex)
    
    
    
    self.front_left_up_explose_vertex=Vertex(self.front_left_up.wx - 2,self.front_left_up.wy - 2,self.front_left_up.wz - 1)
    self.front_left_up_explose_line = self._compute_line_unit(self.front_left_up,self.front_left_up_explose_vertex)
    
    self.back_left_up_explose_vertex=Vertex(self.back_left_up.wx - 2,self.back_left_up.wy - 2,self.back_left_up.wz + 1)
    self.back_left_up_explose_line = self._compute_line_unit(self.back_left_up,self.back_left_up_explose_vertex)
    
    self.back_right_up_explose_vertex=Vertex(self.back_right_up.wx + 2,self.back_right_up.wy - 2,self.back_right_up.wz + 1)
    self.back_right_up_explose_line = self._compute_line_unit(self.back_right_up,self.back_right_up_explose_vertex)
    
    self.front_right_up_explose_vertex=Vertex(self.front_right_up.wx + 2,self.front_right_up.wy - 2,self.front_right_up.wz - 1)
    self.front_right_up_explose_line = self._compute_line_unit(self.front_right_up,self.front_right_up_explose_vertex)
    
    # One cube has 8 corners and every corner generate an explosion line.
    self.explode_lines=[self.front_left_down_explose_line,self.back_left_down_explose_line,self.back_right_down_explose_line,self.front_right_down_explose_line,
                        self.front_left_up_explose_line,self.back_left_up_explose_line,self.back_right_up_explose_line,self.front_right_up_explose_line ]
  
  def explose(self) :
    ''' Cube explosion displaying function. '''
    
    self.explosion_sound_object.play()
    
    max_len=len(self.explode_lines[0])
    idx=0
    size=len(self.explode_lines[0])
    
    
    while idx < max_len :
      # Explosion displaying loop: every corner has an explosion line to follow.
      
      set_point_size(size,False)
      
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
      # Dsipay the playground and the player.
      main_cube.display()
      player.display()
      
      for v in self.explode_lines :
        # Display one item per corner explosion line.
        glColor3ubv(self.self_color)
        glBegin(GL_POINTS)
        glVertex3fv(v[idx].get_vertex())
        glEnd()
        
      sleep(0.15)
      idx += 1
      
      pygame.display.flip()
    
    # Remove the cube instance from the main cube child cubes list.
    main_cube.cubes.remove(self) 
    
    # By an explosion new cube(s) are append to the playground.
    main_cube.generate_cubes()
    
    
  def is_ready_to_explose(self) :
    ''' Check if the cube should explode. '''
    
    if self.life_points == 0 :
      self.update_to_explose()
      self.explose()
    else :
      self.life_points -= 1 
      
    self.is_hit=False  
    sleep(0.15)  
      

class Main_cube_explosing(Base_Class) :
  def _init(self) :
    ''' Class for the final game winning BIG BANG. '''
    
    self.explosion_sound_object=pygame.mixer.Sound("/usr/share/cube-hunter/Sound/Effects/Explosion/Big_bang.wav")
    
    self.size=12
    
    # Explosion centers points:
    self.explosion_start_points=[self.front_left_down,  
				 self.back_left_down,  
				 self.back_right_down,  
				 self.front_right_down, 
				 self.front_left_up,    
				 self.back_left_up,     
				 self.back_right_up,   
				 self.front_right_up,  
                                 self.middle_front_back_at_left_down,  
                                 self.middle_front_back_at_left_up,    
                                 self.center_left_side, 
				 self.middle_front_back_at_right_down,
				 self.middle_front_back_at_right_up,
				 self.center_right_side, 
				 self.middle_left_right_at_front_down, 
				 self.middle_left_right_at_front_up,      
				 self.center_front_side,
				 self.middle_left_right_at_back_down,
				 self.middle_left_right_at_back_up, 
				 self.center_back_side, 
				 self.center_down_side,
				 self.center_up_side,
				 self.cube_center, 
				 self.middle_up_down_at_front_right,
				 self.middle_up_down_at_front_left,
				 self.middle_up_down_at_back_right,
				 self.middle_up_down_at_back_left
				 ]  
    
  def compute_explosion_points(self,start_point) :
    ''' Compute the explosion trajectories. '''
    
    # Every explosion center generate an explosion point in 6 different directions:
    end_point_1=Vertex(start_point.wx,start_point.wy+6,start_point.wz)
    end_point_2=Vertex(start_point.wx,start_point.wy-6,start_point.wz)
    end_point_3=Vertex(start_point.wx,start_point.wy,start_point.wz+6)
    end_point_4=Vertex(start_point.wx,start_point.wy,start_point.wz-6)
    end_point_5=Vertex(start_point.wx+6,start_point.wy,start_point.wz)
    end_point_6=Vertex(start_point.wx-6,start_point.wy,start_point.wz)
    
    res=[]
    res.append(self._compute_line_unit(start_point,end_point_1))
    res.append(self._compute_line_unit(start_point,end_point_2))
    res.append(self._compute_line_unit(start_point,end_point_3))
    res.append(self._compute_line_unit(start_point,end_point_4))
    res.append(self._compute_line_unit(start_point,end_point_5))
    res.append(self._compute_line_unit(start_point,end_point_6))
    
    return res
    
    
  
  def update_to_explose(self) :
    ''' Computing the explosion points directions lines. In relationship to the centers. '''
    
    self.big_bang_lines=[]
    for v in self.explosion_start_points :
      self.big_bang_lines.append(self.compute_explosion_points(v))
    
  def explose(self) :
    ''' Explosion displaying function. '''
    
    # Play the BIG BANG sound.
    self.explosion_sound_object.play()
    
    # Change the displaying speed. 
    sleeper=0.01 
    
    # Set the explosion points color.
    
    
    for v in xrange(0,len(self.big_bang_lines[0][0])) : 
      
      player_win_text=Text(size=1.05,line_width=6,color=(255,0,0))
      player_win_text.construct_string("you win",Vertex(-15.0,-5.0,0.0))
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT) 
      set_point_size(self.size)
      glBegin(GL_POINTS)
      glColor3ubv((127,127,127))
      for y in xrange(0,len(self.big_bang_lines[0])) :
	# 6 points per exploding center displaying from 27 trajectory centers.
	glVertex3fv(self.big_bang_lines[0][y][v].get_vertex())
	glVertex3fv(self.big_bang_lines[1][y][v].get_vertex())
	glVertex3fv(self.big_bang_lines[2][y][v].get_vertex())
	glVertex3fv(self.big_bang_lines[3][y][v].get_vertex())
	glVertex3fv(self.big_bang_lines[4][y][v].get_vertex())
	glVertex3fv(self.big_bang_lines[5][y][v].get_vertex())
	glVertex3fv(self.big_bang_lines[6][y][v].get_vertex())
	glVertex3fv(self.big_bang_lines[7][y][v].get_vertex())
	glVertex3fv(self.big_bang_lines[8][y][v].get_vertex())
	glVertex3fv(self.big_bang_lines[9][y][v].get_vertex())
	glVertex3fv(self.big_bang_lines[10][y][v].get_vertex())
	glVertex3fv(self.big_bang_lines[11][y][v].get_vertex())
	glVertex3fv(self.big_bang_lines[12][y][v].get_vertex())
	glVertex3fv(self.big_bang_lines[13][y][v].get_vertex())
	glVertex3fv(self.big_bang_lines[14][y][v].get_vertex())
	glVertex3fv(self.big_bang_lines[15][y][v].get_vertex())
	glVertex3fv(self.big_bang_lines[16][y][v].get_vertex())
	glVertex3fv(self.big_bang_lines[17][y][v].get_vertex())
	glVertex3fv(self.big_bang_lines[18][y][v].get_vertex())
	glVertex3fv(self.big_bang_lines[19][y][v].get_vertex())
	glVertex3fv(self.big_bang_lines[20][y][v].get_vertex())
	glVertex3fv(self.big_bang_lines[21][y][v].get_vertex())
	glVertex3fv(self.big_bang_lines[22][y][v].get_vertex())
	glVertex3fv(self.big_bang_lines[23][y][v].get_vertex())
	glVertex3fv(self.big_bang_lines[24][y][v].get_vertex())
	glVertex3fv(self.big_bang_lines[25][y][v].get_vertex())
	glVertex3fv(self.big_bang_lines[26][y][v].get_vertex())
      
      glEnd()
      
      glPushMatrix()
      
      camera_view=Localview()
      m=Matrix()                                  # Build the camera position 
      m.translate((0.0,0.0,0.0))                  # Build the camera position 
      m.rotate_vector(0,camera_view.right)        # Build the camera position 
      m.rotate_vector(0,camera_view.up)           # Build the camera position 
      m.rotate_vector(0,camera_view.sight)        # Build the camera position 
      m.load_hardware() 
      
      player_win_text.display()
      
      glPopMatrix()
      
      pygame.display.flip()
      sleep(sleeper)
      
      self.size -= 1
      sleeper += sleeper
    
      
      
    
    
    

class Shoot_ball() :
  def __init__(self,shoot_line,shoot_line_len) :
    ''' Class representing one shoot ball. '''
    
    # Sound to play wennan shoot hit one cube.
    self.shoot_hit_cube_object=pygame.mixer.Sound("/usr/share/cube-hunter/Sound/Effects/shoot_hit_cube/Shoot_hit_cube.wav")
    
    # Shoot line length.
    self.shoot_line_len=shoot_line_len
    
    # Computed shoot ball positions list. 
    self.shoot_line=shoot_line
    
    # Position index in the shoot ball positions list.
    self.idx=0
    
    # Current shoot ball position Vertex.
    self.act_pos=self.shoot_line[0]
    
    # Shoot ball size.
    self.size=16
    
  def display(self) :
    ''' Shoot ball displaying function and test of cube hitting. '''
    
    if self.idx == self.shoot_line_len-1 :
      # The shoot ball end at an intersection.
      try :
        # Removing the shoot ball instance from the player shoot list.
        player.shoot_balls.remove(self)
      except :
	pass
      
      player.is_shooting=False
      
    if not self.idx == self.shoot_line_len : 
      
      ##################################################
      # Display ogf the shoot ball.
      set_point_size(self.size,True)
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT) 
      glBegin(GL_POINTS)
      glVertex3fv(self.shoot_line[self.idx].get_vertex())
      glEnd()
      
      
      #########################################################
      # BUGFIX:
      # This function don't work as wanted because 
      # the shoot ball cube hitting can be not detected or it can be register twice or more times.
      # I have do my best but this bug doesn't not decrement the playing pleasure and the playbility of the game but it appends the bugs mentionned before.
      # and decrement hitten cube life points in relationship of what mentionned before.
      
      for v in main_cube.cubes :
	if not v.is_hit and (self.act_pos.wx < v.act_pos.wx + (1.75 / 2) and self.act_pos.wx > v.act_pos.wx - (1.75 / 2)) and  (self.act_pos.wy < v.act_pos.wy + (1.75 / 2) and self.act_pos.wy > v.act_pos.wy - (1.75 / 2)) and (self.act_pos.wz < v.act_pos.wz + (1.75 / 2) and self.act_pos.wz > v.act_pos.wz - (1.75 / 2)) :
          
          v.is_hit=True  # Try to set an boolean value for not double hitting count.
          
          self.shoot_hit_cube_object.play()
          
          try :
	    # Removing the shoot ball instance from the player shoot list.
            player.shoot_balls.remove(self)
          except :
	    pass
	  
          player.is_shooting=False
          v.is_ready_to_explose()  # Function to check if the cube should explode.
	  return
	
      self.act_pos=self.shoot_line[self.idx]
      
      self.idx += 1

class Player(Base_Class,Player_explosing) :
  def __init__(self,size,color,lifes) :
    
    # Generate base player dislaying datas: triangles, 
    #                                       player middle point computing points list,
    #                                       player cube circonscribe points list for colliding with an cube.
    self.player_display,self.player_center_points,self.collide_center_points=gen_player(1.75)
    #########################################################################################################
    
    self.config_player_explose()  # Inherit exploding settings method from Player_explosing class.
    
    self.self_color=color         # Displaying color data.
    self.size=3                   # Displaying size data. 
    
    self.life_points=5            # life points before exploding.
    self.lifes=lifes              # Can explode so many times before losing. 
    
    self.is_hit=False             # Boolean collide variable.
    self.is_dead=False            # Boolean lose a life variable.
    
    # Control variable for not double down count life points.
    # To BUGFIX because this don't work as wanted.
    self.hit_pos=compute_middle_point(self.collide_center_points)
    
    # Set the player in the right start direction.
    m=Matrix()
    m.rotate_y(90)
    self._update_pos(m)
    
    
    # Inherit method from Base_Class. 
    self.load_shared_datas()  
    
    # Compute slide position vertice.
    self.act_pos=compute_middle_point(self.player_center_points)
    
    # Compute collide with an cube position vertice.
    self.collide_center=compute_middle_point(self.collide_center_points)
    
    m=Matrix()
    res=get_distance_as_vertex(self.act_pos,self.cube_center)
    m.translate((res.wx,res.wy,res.wz))
    self._update_pos(m)
    
    # Set the start player direction.
    self.player_direction='F'
    
        
    
    # Set the start index
    self.pos_index=13
    
    # Set moving and shooting status boolean values.
    self.is_moving=False
    self.is_shooting=False
    
    # List of fire shoots.
    self.shoot_balls=[]
    
        
    
    
    
  
  def load_shared_datas(self) :
    self.initialize_shared_datas()  # Inherit method from Base_Class.
  
  def change_player_direction(self,direction) :
    ''' Function to change the player orientation in relationship to the given direction indicator and the current direction. '''
    
    if self.player_direction == 'F' and direction == 'R' :
      m1=Matrix()
      m1.translate((-self.act_pos.wx,-self.act_pos.wy,-self.act_pos.wz))
      m1.rotate_z(90.)
      m1.translate((self.act_pos.wx,self.act_pos.wy,self.act_pos.wz))
      m2=Matrix()
      center=compute_middle_point(self.collide_center_points)
      res=get_distance_as_vertex(center,self.act_pos)
      m2.translate((res.wx,res.wy,res.wz))
      self._update_pos(m1,m2)
      self.player_direction = 'R'
    
    elif self.player_direction == 'F' and direction == 'L' :
      m1=Matrix()
      m1.translate((-self.act_pos.wx,-self.act_pos.wy,-self.act_pos.wz))
      m1.rotate_z(-90.)
      m1.translate((self.act_pos.wx,self.act_pos.wy,self.act_pos.wz))
      m2=Matrix()
      center=compute_middle_point(self.collide_center_points)
      res=get_distance_as_vertex(center,self.act_pos)
      m2.translate((res.wx,res.wy,res.wz))
      self._update_pos(m1,m2)
      self.player_direction = 'L'
      
    elif self.player_direction == 'F' and direction == 'U' :
      m1=Matrix()
      m1.translate((-self.act_pos.wx,-self.act_pos.wy,-self.act_pos.wz))
      m1.rotate_y(-90.)
      m1.translate((self.act_pos.wx,self.act_pos.wy,self.act_pos.wz))
      m2=Matrix()
      center=compute_middle_point(self.collide_center_points)
      res=get_distance_as_vertex(center,self.act_pos)
      m2.translate((res.wx,res.wy,res.wz))
      self._update_pos(m1,m2)
      self.player_direction = 'U'
      
    elif self.player_direction == 'F' and direction == 'D' :
      m1=Matrix()
      m1.translate((-self.act_pos.wx,-self.act_pos.wy,-self.act_pos.wz))
      m1.rotate_y(90.)
      m1.translate((self.act_pos.wx,self.act_pos.wy,self.act_pos.wz))
      m2=Matrix()
      center=compute_middle_point(self.collide_center_points)
      res=get_distance_as_vertex(center,self.act_pos)
      m2.translate((res.wx,res.wy,res.wz))
      self._update_pos(m1,m2)
      self.player_direction = 'D'
      
      
    elif self.player_direction == 'R' and direction == 'R' :
      m1=Matrix()
      m1.translate((-self.act_pos.wx,-self.act_pos.wy,-self.act_pos.wz))
      m1.rotate_z(90.)
      m1.translate((self.act_pos.wx,self.act_pos.wy,self.act_pos.wz))
      m2=Matrix()
      center=compute_middle_point(self.collide_center_points)
      res=get_distance_as_vertex(center,self.act_pos)
      m2.translate((res.wx,res.wy,res.wz))
      self._update_pos(m1,m2)
      self.player_direction = 'B'
      
    elif self.player_direction == 'R' and direction == 'L' :
      m1=Matrix()
      m1.translate((-self.act_pos.wx,-self.act_pos.wy,-self.act_pos.wz))
      m1.rotate_z(-90.)
      m1.translate((self.act_pos.wx,self.act_pos.wy,self.act_pos.wz))
      m2=Matrix()
      center=compute_middle_point(self.collide_center_points)
      res=get_distance_as_vertex(center,self.act_pos)
      m2.translate((res.wx,res.wy,res.wz))
      self._update_pos(m1,m2)
      self.player_direction = 'F'
      
    elif self.player_direction == 'R' and direction == 'U' :
      m1=Matrix()
      m1.translate((-self.act_pos.wx,-self.act_pos.wy,-self.act_pos.wz))
      m1.rotate_x(90.)
      m1.translate((self.act_pos.wx,self.act_pos.wy,self.act_pos.wz))
      m2=Matrix()
      center=compute_middle_point(self.collide_center_points)
      res=get_distance_as_vertex(center,self.act_pos)
      m2.translate((res.wx,res.wy,res.wz))
      self._update_pos(m1,m2)
      self.player_direction = 'U'
      
    elif self.player_direction == 'R' and direction == 'D' :
      m1=Matrix()
      m1.translate((-self.act_pos.wx,-self.act_pos.wy,-self.act_pos.wz))
      m1.rotate_x(-90.)
      m1.translate((self.act_pos.wx,self.act_pos.wy,self.act_pos.wz))
      m2=Matrix()
      center=compute_middle_point(self.collide_center_points)
      res=get_distance_as_vertex(center,self.act_pos)
      m2.translate((res.wx,res.wy,res.wz))
      self._update_pos(m1,m2)
      self.player_direction = 'D'
      
      
    elif self.player_direction == 'L' and direction == 'R' :
      m1=Matrix()
      m1.translate((-self.act_pos.wx,-self.act_pos.wy,-self.act_pos.wz))
      m1.rotate_z(90.)
      m1.translate((self.act_pos.wx,self.act_pos.wy,self.act_pos.wz))
      m2=Matrix()
      center=compute_middle_point(self.collide_center_points)
      res=get_distance_as_vertex(center,self.act_pos)
      m2.translate((res.wx,res.wy,res.wz))
      self._update_pos(m1,m2)
      self.player_direction = 'F'
      
    elif self.player_direction == 'L' and direction == 'L' :
      m1=Matrix()
      m1.translate((-self.act_pos.wx,-self.act_pos.wy,-self.act_pos.wz))
      m1.rotate_z(-90.)
      m1.translate((self.act_pos.wx,self.act_pos.wy,self.act_pos.wz))
      m2=Matrix()
      center=compute_middle_point(self.collide_center_points)
      res=get_distance_as_vertex(center,self.act_pos)
      m2.translate((res.wx,res.wy,res.wz))
      self._update_pos(m1,m2)
      self.player_direction = 'B'
      
    elif self.player_direction == 'L' and direction == 'U' :
      m=Matrix()
      m.translate((-self.act_pos.wx,-self.act_pos.wy,-self.act_pos.wz))
      m.rotate_x(-90.)
      m.translate((self.act_pos.wx,self.act_pos.wy,self.act_pos.wz))
      self._update_pos(m)
      self.player_direction = 'U'
      
    elif self.player_direction == 'L' and direction == 'D' :
      m1=Matrix()
      m1.translate((-self.act_pos.wx,-self.act_pos.wy,-self.act_pos.wz))
      m1.rotate_x(90.)
      m1.translate((self.act_pos.wx,self.act_pos.wy,self.act_pos.wz))
      m2=Matrix()
      center=compute_middle_point(self.collide_center_points)
      res=get_distance_as_vertex(center,self.act_pos)
      m2.translate((res.wx,res.wy,res.wz))
      self._update_pos(m1,m2)
      self.player_direction = 'D'
      
    
    elif self.player_direction == 'B' and direction == 'R' :
      m1=Matrix()
      m1.translate((-self.act_pos.wx,-self.act_pos.wy,-self.act_pos.wz))
      m1.rotate_z(90.)
      m1.translate((self.act_pos.wx,self.act_pos.wy,self.act_pos.wz))
      m2=Matrix()
      center=compute_middle_point(self.collide_center_points)
      res=get_distance_as_vertex(center,self.act_pos)
      m2.translate((res.wx,res.wy,res.wz))
      self._update_pos(m1,m2)
      self.player_direction = 'L'
      
    elif self.player_direction == 'B' and direction == 'L' :
      m1=Matrix()
      m1.translate((-self.act_pos.wx,-self.act_pos.wy,-self.act_pos.wz))
      m1.rotate_z(-90.)
      m1.translate((self.act_pos.wx,self.act_pos.wy,self.act_pos.wz))
      m2=Matrix()
      center=compute_middle_point(self.collide_center_points)
      res=get_distance_as_vertex(center,self.act_pos)
      m2.translate((res.wx,res.wy,res.wz))
      self._update_pos(m1,m2)
      self.player_direction = 'R'
      
    elif self.player_direction == 'B' and direction == 'U' :
      m1=Matrix()
      m1.translate((-self.act_pos.wx,-self.act_pos.wy,-self.act_pos.wz))
      m1.rotate_y(90.)
      m1.translate((self.act_pos.wx,self.act_pos.wy,self.act_pos.wz))
      m2=Matrix()
      center=compute_middle_point(self.collide_center_points)
      res=get_distance_as_vertex(center,self.act_pos)
      m2.translate((res.wx,res.wy,res.wz))
      self._update_pos(m1,m2)
      self.player_direction = 'U'
      
    elif self.player_direction == 'B' and direction == 'D' :
      m1=Matrix()
      m1.translate((-self.act_pos.wx,-self.act_pos.wy,-self.act_pos.wz))
      m1.rotate_y(-90.)
      m1.translate((self.act_pos.wx,self.act_pos.wy,self.act_pos.wz))
      m2=Matrix()
      center=compute_middle_point(self.collide_center_points)
      res=get_distance_as_vertex(center,self.act_pos)
      m2.translate((res.wx,res.wy,res.wz))
      self._update_pos(m1,m2)
      self.player_direction = 'D'
       
    
    elif self.player_direction == 'U' and direction == 'U' :
      m1=Matrix()
      m1.translate((-self.act_pos.wx,-self.act_pos.wy,-self.act_pos.wz))
      m1.rotate_y(90.)
      m1.translate((self.act_pos.wx,self.act_pos.wy,self.act_pos.wz))
      m2=Matrix()
      center=compute_middle_point(self.collide_center_points)
      res=get_distance_as_vertex(center,self.act_pos)
      m2.translate((res.wx,res.wy,res.wz))
      self._update_pos(m1,m2)
      self.player_direction = 'F'
      
    elif self.player_direction == 'U' and direction == 'D' :
      m1=Matrix()
      m1.translate((-self.act_pos.wx,-self.act_pos.wy,-self.act_pos.wz))
      m1.rotate_y(-90.)
      m1.translate((self.act_pos.wx,self.act_pos.wy,self.act_pos.wz))
      m2=Matrix()
      center=compute_middle_point(self.collide_center_points)
      res=get_distance_as_vertex(center,self.act_pos)
      m2.translate((res.wx,res.wy,res.wz))
      self._update_pos(m1,m2)
      self.player_direction = 'B'
      
    elif self.player_direction == 'U' and direction == 'L' :
      m1=Matrix()
      m1.translate((-self.act_pos.wx,-self.act_pos.wy,-self.act_pos.wz))
      m1.rotate_x(90.)
      m1.translate((self.act_pos.wx,self.act_pos.wy,self.act_pos.wz))
      m2=Matrix()
      center=compute_middle_point(self.collide_center_points)
      res=get_distance_as_vertex(center,self.act_pos)
      m2.translate((res.wx,res.wy,res.wz))
      self._update_pos(m1,m2)
      self.player_direction = 'L'
      
      
    elif self.player_direction == 'U' and direction == 'R' :
      m1=Matrix()
      m1.translate((-self.act_pos.wx,-self.act_pos.wy,-self.act_pos.wz))
      m1.rotate_x(-90.)
      m1.translate((self.act_pos.wx,self.act_pos.wy,self.act_pos.wz))
      m2=Matrix()
      center=compute_middle_point(self.collide_center_points)
      res=get_distance_as_vertex(center,self.act_pos)
      m2.translate((res.wx,res.wy,res.wz))
      self._update_pos(m1,m2)
      self.player_direction = 'R'
      
    
    
    elif self.player_direction == 'D' and direction == 'U' :
      m1=Matrix()
      m1.translate((-self.act_pos.wx,-self.act_pos.wy,-self.act_pos.wz))
      m1.rotate_y(-90.)
      m1.translate((self.act_pos.wx,self.act_pos.wy,self.act_pos.wz))
      m2=Matrix()
      center=compute_middle_point(self.collide_center_points)
      res=get_distance_as_vertex(center,self.act_pos)
      m2.translate((res.wx,res.wy,res.wz))
      self._update_pos(m1,m2)
      self.player_direction = 'F'
      
    elif self.player_direction == 'D' and direction == 'D' :
      m1=Matrix()
      m1.translate((-self.act_pos.wx,-self.act_pos.wy,-self.act_pos.wz))
      m1.rotate_y(90.)
      m1.translate((self.act_pos.wx,self.act_pos.wy,self.act_pos.wz))
      m2=Matrix()
      center=compute_middle_point(self.collide_center_points)
      res=get_distance_as_vertex(center,self.act_pos)
      m2.translate((res.wx,res.wy,res.wz))
      self._update_pos(m1,m2)
      self.player_direction = 'B'
      
    elif self.player_direction == 'D' and direction == 'L' :
      m1=Matrix()
      m1.translate((-self.act_pos.wx,-self.act_pos.wy,-self.act_pos.wz))
      m1.rotate_x(-90.)
      m1.translate((self.act_pos.wx,self.act_pos.wy,self.act_pos.wz))
      m2=Matrix()
      center=compute_middle_point(self.collide_center_points)
      res=get_distance_as_vertex(center,self.act_pos)
      m2.translate((res.wx,res.wy,res.wz))
      self._update_pos(m1,m2)
      self.player_direction = 'L'
      
      
    elif self.player_direction == 'D' and direction == 'R' :
      m1=Matrix()
      m1.translate((-self.act_pos.wx,-self.act_pos.wy,-self.act_pos.wz))
      m1.rotate_x(90.)
      m1.translate((self.act_pos.wx,self.act_pos.wy,self.act_pos.wz))
      m2=Matrix()
      center=compute_middle_point(self.collide_center_points)
      res=get_distance_as_vertex(center,self.act_pos)
      m2.translate((res.wx,res.wy,res.wz))
      self._update_pos(m1,m2)
      self.player_direction = 'R'
      
  def lose_a_life(self) :
    # Reset the variable by losing one life.
    self.player_display=[]
    self.is_hit=True
   
  def new_life(self) :
    # Reset the settings play a new life.
    self.is_dead=False
    self.is_hit=False
    self.player_display,self.player_center_points,self.collide_center_points=gen_player(1.75)
    self.hit_pos=compute_middle_point(self.collide_center_points)
    
    # Set the player in the right start direction.
    m=Matrix()
    m.rotate_y(90)
    self._update_pos(m)
    
    # Compute slide position vertice.
    self.act_pos=compute_middle_point(self.player_center_points)
    
    # Compute collide with an cube position vertice.
    self.collide_center=compute_middle_point(self.collide_center_points)
    
    # Set the start player direction.
    self.player_direction='F'
    
    # Set the start index
    self.pos_index=13
    
    # Set moving and shooting status boolean values.
    self.is_moving=False
    self.is_shooting=False
    
    # List of fire shoots.
    self.shoot_balls=[]
      
  def collide_check(self) :
    ''' Function to check an collision with an cube. '''
    
    #########################################################
    # BUGFIX:
    # This function don't work as wanted because 
    # the colliding can be not detected or it can be register twice or more times.
    # I have do my best but this bug doesn't not decrement the playing pleasure and the playbility of the game it only false the lifes points count.
     
    for v in main_cube.cubes :
      
      if not self.is_hit and ( v.center.wx > self.collide_center.wx - (1.75/2) and v.center.wy > self.collide_center.wy - (1.75/2) and v.center.wz > self.collide_center.wz - (1.75/2) ) and ( v.center.wx < self.collide_center.wx + (1.75/2) and v.center.wy < self.collide_center.wy + (1.75/2) and v.center.wz < self.collide_center.wz + (1.75/2) ) :
	if self.hit_pos.get_vertex() == self.collide_center.get_vertex() :
	  # Temporar colliding position check.
	  return
	
	self.is_hit=True                  # An variable for the bugfix try.
	self.hit_pos=self.collide_center  # An variable for the bugfix try. 
	
	player.is_ready_to_explose()      # Check if the player lose a life and explode.
  
  def shoot(self) :
    ''' Main fire function '''
    if not self.is_moving and not self.is_shooting and not self.is_dead :
      
      self.is_shooting = True
      
      # Start shoot ball position.
      shoot_from  = self.act_pos
      
      # End shoot ball position.
      shoot_end   = self.position[self.pos_index].get(self.player_direction)
      
      if shoot_end :
	
        if self.position[shoot_end[1]].get(self.player_direction) :
          # Get the shoot destination position if it doesn't hit an cube.   
          shoot_end = self.position[shoot_end[1]].get(self.player_direction)[0]
        else :
	  # Get the shoot destination position if it doesn't hit an cube.  
          shoot_end   = self.position[self.pos_index].get(self.player_direction)[0]
        
        # Compute the shoot line.
        shoot_line_len, shoot_line = self._compute_shoot_line(shoot_from,shoot_end)
        
        # Append the shoot ball to container.
        self.shoot_balls.append(Shoot_ball(shoot_line,shoot_line_len))
        
      else :
        self.is_shooting=False
        
  def _compute_shoot_line(self,vertex1,vertex2) :
    ''' Function to compute every position from the shoot line. '''
    
    res=[]
    res.append(vertex1)
    res.append(get_middle(vertex1,vertex2))
    res.append(vertex2)
    
    
    limit= int(get_distance(vertex1,vertex2))  
    breaker=0

    while breaker < limit :
      breaker,res=self._computing_line_unit(res)
      
    return breaker,res
    
  
  def _computing_shoot_line_unit(self,res) :
    ''' Function to compute one position from the shoot line and the break of the computing. '''
    tmp_res=res
    res=res[0::]
    i=0
    ii=1
    while i < len(res) :
      if not i == len(res)-1 :
        tmp_res.insert(i+ii,get_middle(res[i],res[i+1]))
        
      else :
        break
      i += 1
      ii += 1
      
    res=tmp_res
    
    return len(tmp_res),res   
  
  def move(self) :
    ''' Main move forward function. '''
    
    if self.is_moving :
      # The player is moving.
      return
    
    if self.position[self.pos_index].get(self.player_direction) :
      
      # Start position from the player slide: at an intersection.
      self.from_position = self.position[self.pos_index].get('pos')[0]
      
      # End position from the player slide taken from the dict inherit from Base_Class.
      self.to_position   = self.position[self.pos_index].get(self.player_direction)[0]
      
      # Index from the player end slide position in the list of moving dicts contains in the array inherit from Base_Class.
      self.pos_index     = self.position[self.pos_index].get(self.player_direction)[1]
      
      # Player slide steps computing.
      moveto=self._compute_line_unit(self.from_position,self.to_position)
      
      self.is_moving=True
      
      self.slide_to(moveto)
    
    
  
  def slide_to(self,slide_to) :
    if slide_to :
      # Set sliding values.
      self.slide_line     = slide_to
      self.slide_line_idx = 0
      self.slide_line_len = len(slide_to)
        
  def slide(self) :
    ''' Main plyetr sliding function. '''
    
    if self.is_moving and not self.is_shooting :
      
      if self.slide_line_idx + 1 < self.slide_line_len :
        self.slide_line_idx  += 1  # Increment steps index in the sliding positions list.
        self.act_pos=self.slide_line[self.slide_line_idx] # Set the current position from the slide list.
        
        
        
        # Update player slide position.
        m1=Matrix()
        center=compute_middle_point(self.player_center_points)
        res=get_distance_as_vertex(center,self.act_pos)
        m1.translate((res.wx,res.wy,res.wz))
        
        # Update collide point. 
        m2=Matrix()
        center=compute_middle_point(self.collide_center_points)
        res=get_distance_as_vertex(center,self.act_pos)
        m2.translate((res.wx,res.wy,res.wz))
        
        # Effective player settings updating.
        self._update_pos(m1,m2)
        
      else :
        self.is_moving=False
    
  def display(self) :
    ''' Player display function. '''
    
    glColor3ubv(self.self_color)
    set_lined(self.size)
    
    for v in self.player_display :
      glBegin(GL_TRIANGLES)
      for y in v :
	glVertex3fv(y.get_vertex())
      glEnd()
  
  def _update_pos(self,matrix1,matrix2=False) :
    ''' Player slide and colliding positions updates. '''
    
    if not self.is_dead :
      # Player slide position update. 
      
      tmp=[]
      for v in self.player_center_points :
	tmp.append(matrix1.mult_vertex(v))
	
      self.player_center_points=tmp
      self.act_pos = compute_middle_point(tmp)
      
      if matrix2 :
	# Player collide position update. 
	tmp=[]
	for v in self.collide_center_points :
	  tmp.append(matrix2.mult_vertex(v))
	
	self.collide_center_points=tmp
	self.collide_center = compute_middle_point(tmp)
      
      self._update_coords()
	
      
    
    
    
    
  
  def _update_coords(self) :
    ''' Player display datas updating. '''
    
    pointed=self.player_center_points.pop(0)
    
    self.player_display=[]
    i=0
    while i < len(self.player_center_points) :
      # Loop to build triangles for player pentahedron displaying.
      tmp=[]
      if i+1 != len(self.player_center_points) : 
        tmp.append(self.player_center_points[i])
        tmp.append(self.player_center_points[i+1])
        tmp.append(pointed)
      else :
        tmp.append(self.player_center_points[i])
        tmp.append(self.player_center_points[0])
        tmp.append(pointed)  
      self.player_display.append(tmp)   
      i += 1 
    
    self.player_center_points.insert(0,pointed)
    
    
    

    
    





class Cube_IA(Base_Class,Cube_explosing) :
  
  def __init__(self,color,life_points,game_on) :
    # Display base settings.
    side_length=1.75
    self.self_vertex_list=[Vertex(-side_length/2.,-side_length/2., -side_length/2.), Vertex(side_length/2.,-side_length/2., -side_length/2.), Vertex(side_length/2., side_length/2., -side_length/2.), Vertex(-side_length/2., side_length/2., -side_length/2.),
                           Vertex(-side_length/2.,-side_length/2.,  side_length/2.), Vertex(side_length/2.,-side_length/2.,  side_length/2.), Vertex(side_length/2., side_length/2.,  side_length/2.), Vertex(-side_length/2., side_length/2.,  side_length/2.)]
    
    
    self.life_points=life_points  # Number of life points before exploding.
    self.is_hit=False             # Boolean value case the cube is hitten by an shoot from the player.
    
    if life_points >= 3 and life_points <= 5 :
      self.mode="filled"
      self.self_border_color=(randint(127,255),randint(127,255),randint(127,255))
    else :
      self.mode="lined"
    
    self.size=3
    self.self_color=color
    
    self.initialize_shared_datas()    # Inherit method from Base_Class.
    self.config_cube_explose()        # Inherit method from Cube_explosing class. 
    
    self._update_coords()             # Update the datas needed for cube displaying
    
      
    
    
    
    
    # Random start position.
    random_start_position_from=choice(self.position)  
    
    if game_on :
    
      global player
      
      while random_start_position_from.get('pos')[0].get_vertex() == player.collide_center.get_vertex() :
	random_start_position_from=choice(self.position) 
     
     
    
    # Random destination position from the Base_Class self.start_direction_destination_random variable.
    random_start_position_to = self.position[choice(self.start_direction_destination_random.get(self.position.index(random_start_position_from)))] 
    
    self.slide_to_idx=random_start_position_to.get('pos')[1]    # Index in the look up from the self.position containing possible moves dicts from the Base_Class. 
    
    # Compute the first trajectory.
    self.slide_line=self._compute_line_unit(random_start_position_from.get('pos')[0],random_start_position_to.get('pos')[0])   
    
    
    self.center = self.slide_line[0]  # Vertice of actual postion needed
    
    
    
    self.slide_line_len = len(self.slide_line)-1
    self.slide_line_idx = 0
    
    self.act_pos=self.slide_line[0]        # Vertice needed for next step computing.
    
    m=Matrix()
    m.translate((self.act_pos.wx,self.act_pos.wy,self.act_pos.wz))
    self._update_pos(m)
    
      
      
  
        
  
  def set_new_position_to_go(self) :    
    new_dest=False
    while not new_dest :  # Because we can get an impossible direction to go.
      
      dest_orientation=choice(['F','B','R','L','U','D']) # Take an random direction. 
      
      if self.position[self.slide_to_idx].get(dest_orientation) :
        
        new_dest=self.position[self.slide_to_idx].get(dest_orientation) # New trajectory destination. 
        
        self.slide_line=self._compute_line_unit(self.position[self.slide_to_idx].get('pos')[0],new_dest[0])
        
        self.slide_line_idx=0         # Set the index of current trajectory on start position.
        
        self.slide_to_idx=new_dest[1] # Set next destination index in the look up dict contains in Base_Class self.position variable.
        
  def slide(self) :
    ''' cube slide function: compute and update new step if not new random trajectory. '''
    
    try :
      # By index error it fails, it's secure than an if statement. 
      self.act_pos=self.slide_line[self.slide_line_idx]
    except :
      # Compute new trajectory.
      self.set_new_position_to_go() 
      self.act_pos=self.slide_line[self.slide_line_idx]
    
    # Compute next step to go.
    m=Matrix()
    res=get_distance_as_vertex(self.center,self.act_pos)
    m.translate((res.wx,res.wy,res.wz))
    self._update_pos(m)
    
    self.display()
    
    self.slide_line_idx  += 1
  

  def _update_pos(self,matrix) :
    ''' updating the displaying datas. '''
    tmp_1=[]
    
    for v in self.self_vertex_list :
      tmp_1.append(matrix.mult_vertex(v))
    
    
    self.self_vertex_list=tmp_1
    self.center=compute_middle_point(tmp_1)
    
    self._update_coords()
    
  def _update_coords(self) :
    ''' Displaying datas updating function. ''' 
    self.front = [self.self_vertex_list[0],self.self_vertex_list[1],self.self_vertex_list[2],self.self_vertex_list[3]]
    self.back  = [self.self_vertex_list[4],self.self_vertex_list[5],self.self_vertex_list[6],self.self_vertex_list[7]]
    
    
    self.left  = [self.self_vertex_list[7],self.self_vertex_list[4],self.self_vertex_list[0],self.self_vertex_list[3]]
    self.right = [self.self_vertex_list[5],self.self_vertex_list[6],self.self_vertex_list[1],self.self_vertex_list[2]]
    
    self.up    = [self.self_vertex_list[4],self.self_vertex_list[5],self.self_vertex_list[1],self.self_vertex_list[0]]
    self.down  = [self.self_vertex_list[7],self.self_vertex_list[6],self.self_vertex_list[2],self.self_vertex_list[3]]
      
    
    self.quads_lined_base   = [self.front,self.back]
    self.quads_filled_base  = [self.front,self.back,self.left,self.right,self.up,self.down]
    
  def display(self) :
    
    glColor3ubv(self.self_color)
    
    if self.mode == 'lined' :
      set_lined(self.size)
      
      for v in self.quads_lined_base :
	# We display 2 quads.
        glBegin(GL_LINE_LOOP)
        for y in v :
          glVertex3fv(y.get_vertex())
        glEnd()  
      
      
      for y in range(0,4) :
	# We join the 2 quads.
        glBegin(GL_LINES)
        glVertex3fv(self.quads_lined_base[0][y].get_vertex())
        glVertex3fv(self.quads_lined_base[1][y].get_vertex())
        glEnd()
      
            
          
          
    if self.mode == "filled" :
      set_polygon_mode(mode=GL_FILL)
      
      for v in self.quads_filled_base :
	# Displaying every cube face.
        glBegin(GL_QUADS) 
        for y in v :
          glVertex3fv(y.get_vertex())
        glEnd() 
      
      set_lined(self.size)
      glColor3ubv(self.self_border_color)
      for v in self.quads_lined_base :
	# Displaying every cube border.
        glBegin(GL_LINE_LOOP)
        for y in v :
          glVertex3fv(y.get_vertex())
        glEnd()  
      
      
      for y in range(0,4) :
        glBegin(GL_LINES)
        glVertex3fv(self.quads_lined_base[0][y].get_vertex())
        glVertex3fv(self.quads_lined_base[1][y].get_vertex())
        glEnd()  
    



class Main_Cube(Main_Cube_Base_Class) :
  
  def __init__(self,color,max_childs=5,total_childs=25,size=3) :
     
    self.self_color=color 
    
    self.size=size
    
    self.number_of_childs_to_add=randint(1,max_childs/2+1)
    self.max_childs=max_childs
    self.total_childs=total_childs
    self.death_child_counter=0
    
    self.cubes_counter=0
    
    self.cubes=[]
    
    self.load_shared_datas() 
    
    
  def load_shared_datas(self) :
    self.initialize_shared_datas()  # Method inherit from the base class. 
    ############################################################################################################################################################################
    # Sequences of Quads:
    #
    # Front-Back: 
    self.quad_middles_front_back=[self.middle_front_back_at_left_down,self.middle_front_back_at_right_down,self.middle_front_back_at_right_up,self.middle_front_back_at_left_up] 
    # Left-Right:
    self.quad_middles_left_right=[self.middle_left_right_at_front_down,self.middle_left_right_at_back_down,self.middle_left_right_at_back_up,self.middle_left_right_at_front_up] 
    # Up-Down: 
    self.quad_middles_up_down=[self.middle_up_down_at_front_right,self.middle_up_down_at_front_left,self.middle_up_down_at_back_left,self.middle_up_down_at_back_right] 
    ############################################################################################################################################################################
    
    
    ############################################################################################################################################################################
    # Sequences of Lines:
    #
    # Right-Left centers:
    self.centers_left_right = [self.center_left_side,self.center_right_side,self.center_left_side,self.center_right_side]
    # Front-Back:
    self.centers_front_back = [self.center_front_side,self.center_back_side,self.center_front_side,self.center_back_side]
    # Up-Down:
    self.centers_up_down    = [self.center_down_side,self.center_up_side,self.center_down_side,self.center_up_side]
    ############################################################################################################################################################################
    
    ######################################################################################################################################################################################
    # Drawing containers:                                                                                                                                                             
    #                                                                                                                                                                                 
    self.quads        = [(self.front,self.back),(self.center_left_side,self.center_right_side),(self.center_front_side,self.center_back_side),(self.center_down_side,self.center_up_side)]
    self.cross_quads  = [self.quad_middles_front_back,self.quad_middles_left_right,self.quad_middles_up_down]                                                                         
    ######################################################################################################################################################################################

  def init_inside_room(self) :
    ''' Childs cubes generating at start. '''
    
    for v in xrange(0,self.number_of_childs_to_add) :
      self.add_ia_cube((randint(127,255),randint(127,255),randint(127,255)),randint(1,6))
      self.cubes_counter += 1
      
    self._config_number_of_child_cubes()  
    
    return
    
    
  def generate_cubes(self) :
    ''' Generating main cube child(s) cubes after one child has explosed. '''
    
    self.death_child_counter += 1
    
    if not self.death_child_counter >= self.total_childs :  # The child cubes number is not reached totals.
      
      if self.number_of_childs_to_add != 0 and not len(self.cubes) == 0 :
	# Case we add childs cubes and the playground contains some child cubes. 
	
	for v in xrange(0,self.number_of_childs_to_add) :
	  self.add_ia_cube((randint(127,255),randint(127,255),randint(127,255)),randint(1,6))
      
      elif self.number_of_childs_to_add == 0 and len(self.cubes) == 0 :
	# Case we cannot add some child cubes and the playground contains no child cubes.
	# We must add one because the playground must contains some child cubes so that by
	# an explosion we can generate others.
	
	for v in xrange(0,1) :
	  self.add_ia_cube((randint(127,255),randint(127,255),randint(127,255)),randint(1,6))
      
      elif self.number_of_childs_to_add != 0 and len(self.cubes) == 0 :
        # Case we can add childs without problems.	
	
	for v in xrange(0,self.number_of_childs_to_add) :
	  self.add_ia_cube((randint(127,255),randint(127,255),randint(127,255)),randint(1,6))
	
      # The else statement in case we don't add some child is not needed. 	
      main_cube._config_number_of_child_cubes() 
      
	
    else :
      
      if len(self.cubes) == 0 :
	final_explosion=Main_cube_explosing()
	final_explosion.initialize_shared_datas()
	final_explosion._init()
	final_explosion.update_to_explose()
	final_explosion.explose()
	event_to_send=pygame.event.Event(USEREVENT,{"code":"end_game"})
	pygame.event.post(event_to_send)
	
    
    
  def _config_number_of_child_cubes(self) :
    ''' Set the number of main cubes child cubes to add at next explosion. '''
    to_generate=randint(0,abs(self.max_childs-len(self.cubes))+1)
    if to_generate > (self.total_childs-self.cubes_counter) :
      to_generate=to_generate-self.cubes_counter
      if to_generate < 0 :
	to_generate=0
	
    self.number_of_childs_to_add=to_generate  
    self.cubes_counter += self.number_of_childs_to_add 
      
  
  def add_ia_cube(self,color,life_points) :
    ''' Function to random the start position and trajectory from an added child cube. '''
    
    global player
    
    
    
    # Instantiate an child cube with random settings.
    cube_ia=Cube_IA(color,life_points,True) 
       
    # Append the new child cube to the main cubes childs list.
    self.cubes.append(cube_ia)
    
  def manage_inside_room(self) :
    ''' Main child cubes moving function. '''
    for v in self.cubes :
      v.slide()  
    
  
  def display(self) :
    ''' Main cube display function. '''
    
    glColor3ubv(self.self_color)

    set_lined(self.size)
      
      
    for v in self.quads[0] :
      # Display 2 big quads.
      glBegin(GL_POLYGON)
      for y in v :
	glVertex3fv(y.get_vertex())
      glEnd()  
    
    
    for y in range(0,4) :
      # Join the 2 big quads.
      glBegin(GL_LINES)
      glVertex3fv(self.quads[0][0][y].get_vertex())
      glVertex3fv(self.quads[0][1][y].get_vertex())
      glEnd()
    
    
    
    for v in self.cross_quads :
      # Display the cross from every cube side.
      glBegin(GL_QUADS)
      for y in v :
	glVertex3fv(y.get_vertex())
      glEnd()
    
    i=-1
    while i > -4 :
      # Join the the centers from every side.
      glBegin(GL_LINES)
      for v in self.quads[i] :
	glVertex3fv(v.get_vertex())
      glEnd()
      i -= 1

class Text(object) :
  def __init__(self,size=0.75,line_width=3,color=(255,0,0)) :
    
    self.line_width=line_width
    self.color=color
    self.size=size
    
    self.vertex_list_a=[Vertex(-size,size*3.0,0.0),Vertex(-size,-size,0.0),Vertex(size,-size,0.0),Vertex(size,size*3.0,0.0),Vertex(size,size,0.0),Vertex(-size,size,0.0)]
    self.vertex_list_b=[Vertex(-size,-size,0.0),Vertex(-size,size*3.0,0.0),Vertex(size,size*3.0,0.0),Vertex(size,size,0.0),Vertex(-size,size,0.0)]
    self.vertex_list_c=[Vertex(size,-size,0.0),Vertex(-size,-size,0.0),Vertex(-size,size*3.0,0.0),Vertex(size,size*3.0,0.0)]
    self.vertex_list_d=[Vertex(-size,-size,0.0),Vertex(-size,size*3.0,0.0),Vertex(0.0,size*3.0,0.0),Vertex(size,size*1.5,0.0),Vertex(size,size*0.5,0.0),Vertex(0.0,-size,0.0),Vertex(-size,-size,0.0)]
    self.vertex_list_e=[Vertex(size,size*3.0,0.0),Vertex(-size,size*3.0,0.0),Vertex(-size,size,0.0),Vertex(size,size,0.0),Vertex(-size,size,0.0),Vertex(-size,-size,0.0),Vertex(size,-size,0.0)]
    self.vertex_list_f=[Vertex(-size,size*3.0,0.0),Vertex(-size,size,0.0),Vertex(size,size,0.0),Vertex(-size,size,0.0),Vertex(-size,-size,0.0),Vertex(size,-size,0.0)]
    self.vertex_list_g=[Vertex(size,-size,0.0),Vertex(-size,-size,0.0),Vertex(-size,size*3.0,0.0),Vertex(size,size*3.0,0.0),Vertex(size,size,0.0),Vertex(0.0,size,0.0)]
    self.vertex_list_h=[Vertex(-size,-size,0.0),Vertex(-size,size*3.0,0.0),Vertex(-size,size,0.0),Vertex(size,size,0.0),Vertex(size,size*3.0,0.0),Vertex(size,-size,0.0)]
    self.vertex_list_i=[Vertex(-size,-size,0.0),Vertex(size,-size,0.0),Vertex(0.0,-size,0.0),Vertex(0.0,size*3.0,0.0),Vertex(-size,size*3.0,0.0),Vertex(size,size*3.0,0.0)]
    self.vertex_list_j=[Vertex(-size,-size,0.0),Vertex(size,-size,0.0),Vertex(0.0,-size,0.0),Vertex(0.0,size*3.0,0.0),Vertex(-size,size*3.0,0.0)] 
    self.vertex_list_k=[Vertex(-size,-size,0.0),Vertex(-size,size*3.0,0.0),Vertex(-size,size,0.0),Vertex(size,-size,0.0),Vertex(-size,size,0.0),Vertex(size,size*3.0,0.0)] 
    self.vertex_list_l=[Vertex(-size,-size,0.0),Vertex(-size,size*3.0,0.0),Vertex(size,size*3.0,0.0)]
    self.vertex_list_m=[Vertex(-size,size*3.0,0.0),Vertex(-size,-size,0.0),Vertex(0.0,size,0.0),Vertex(size,-size,0.0),Vertex(size,size*3.0,0.0)] 
    self.vertex_list_n=[Vertex(-size,size*3.0,0.0),Vertex(-size,-size,0.0),Vertex(size,size*3.0,0.0),Vertex(size,-size,0.0)] 
    self.vertex_list_o=[Vertex(-size,size*3.0,0.0),Vertex(-size,-size,0.0),Vertex(size,-size,0.0),Vertex(size,size*3.0,0.0),Vertex(-size,size*3.0,0.0)] 
    self.vertex_list_p=[Vertex(-size,size*3.0,0.0),Vertex(-size,-size,0.0),Vertex(size,-size,0.0),Vertex(size,size,0.0),Vertex(-size,size,0.0)] 
    self.vertex_list_q=[Vertex(size,size*3.0,0.0),Vertex(size,-size,0.0),Vertex(-size,-size,0.0),Vertex(-size,size,0.0),Vertex(size,size,0.0)]
    self.vertex_list_r=[Vertex(-size,size*3.0,0.0),Vertex(-size,-size,0.0),Vertex(size,-size,0.0)] 
    self.vertex_list_s=[Vertex(-size,size*3.0,0.0),Vertex(size,size*3.0,0.0),Vertex(size,size,0.0),Vertex(-size,size,0.0),Vertex(-size,-size,0.0),Vertex(size,-size,0.0)] 
    self.vertex_list_t=[Vertex(-size,-size,0.0),Vertex(size,-size,0.0),Vertex(0.0,-size,0.0),Vertex(0.0,size*3.0,0.0)] 
    self.vertex_list_u=[Vertex(-size,-size,0.0),Vertex(-size,size*3.0,0.0),Vertex(size,size*3.0,0.0),Vertex(size,-size,0.0)] 
    self.vertex_list_v=[Vertex(-size,-size,0.0),Vertex(0.0,size*3.0,0.0),Vertex(size,-size,0.0)] 
    self.vertex_list_w=[Vertex(-size,-size,0.0),Vertex(-size,size*3.0,0.0),Vertex(0.0,size,0.0),Vertex(size,size*3.0,0.0),Vertex(size,-size,0.0)]
    self.vertex_list_x=[Vertex(-size,-size,0.0),Vertex(size,size*3.0,0.0),Vertex(0.0,size,0.0), Vertex(size,-size,0.0), Vertex(-size,size*3.0,0.0)] 
    self.vertex_list_y=[Vertex(0.0,size*3.0,0.0),Vertex(0.0,size,0.0),Vertex(size,-size,0.0),Vertex(0.0,size,0.0), Vertex(-size,-size,0.0)] 
    self.vertex_list_z=[Vertex(-size,-size,0.0),Vertex(size,-size,0.0),Vertex(-size,size*3.0,0.0), Vertex(size,size*3.0,0.0)] 
    
    self.vertex_list_0=[Vertex(-size,size*3.0,0.0),Vertex(-size,-size,0.0),Vertex(size,-size,0.0),Vertex(size,size*3.0,0.0),Vertex(-size,size*3.0,0.0)] 
    self.vertex_list_1=[Vertex(-size,-size,0.0),Vertex(size,-size,0.0),Vertex(size,size*3.0,0.0)]
    self.vertex_list_2=[Vertex(size,size*3.0,0.0),Vertex(-size,size*3.0,0.0),Vertex(-size,size,0.0),Vertex(size,size,0.0),Vertex(size,-size,0.0),Vertex(-size,-size,0.0)]
    self.vertex_list_3=[Vertex(-size,size*3.0,0.0),Vertex(size,size*3.0,0.0),Vertex(size,size,0.0),Vertex(-size,size,0.0),Vertex(size,size,0.0),Vertex(size,-size,0.0),Vertex(-size,-size,0.0)]
    self.vertex_list_4=[Vertex(size,size*3.0,0.0),Vertex(size,-size,0.0),Vertex(size,size,0.0), Vertex(-size,size,0.0), Vertex(-size,-size,0.0)]
    self.vertex_list_5=[Vertex(-size,size*3.0,0.0),Vertex(size,size*3.0,0.0),Vertex(size,size,0.0),Vertex(-size,size,0.0),Vertex(-size,-size,0.0),Vertex(size,-size,0.0)] 
    self.vertex_list_6=[Vertex(size,-size,0.0),Vertex(-size,-size,0.0),Vertex(-size,size*3.0,0.0),Vertex(size,size*3.0,0.0),Vertex(size,size,0.0),Vertex(-size,size,0.0)]
    self.vertex_list_7=[Vertex(size,size*3.0,0.0),Vertex(size,size,0.0),Vertex(-size,size,0.0),Vertex(size,size,0.0),Vertex(size,-size,0.0),Vertex(-size,-size,0.0)]
    self.vertex_list_8=[Vertex(-size,size*3.0,0.0),Vertex(-size,-size,0.0),Vertex(size,-size,0.0),Vertex(size,size*3.0,0.0),Vertex(-size,size*3.0,0.0),Vertex(-size,size,0.0),Vertex(size,size,0.0)] 
    self.vertex_list_9=[Vertex(-size,size*3.0,0.0),Vertex(size,size*3.0,0.0),Vertex(size,-size,0.0),Vertex(-size,-size,0.0),Vertex(-size,size,0.0),Vertex(size,size,0.0)]
    
    self.vertex_list_colons=[Vertex(0.0,size*0.25,0.0),Vertex(0.0,size*2.25,0.0)]
    
    
    self.string_look_up={'a': self.vertex_list_a,
			 'b': self.vertex_list_b,
			 'c': self.vertex_list_c,
			 'd': self.vertex_list_d,
			 'e': self.vertex_list_e,
			 'f': self.vertex_list_f,
			 'g': self.vertex_list_g,
			 'h': self.vertex_list_h,
			 'i': self.vertex_list_i,
			 'j': self.vertex_list_j,
			 'k': self.vertex_list_k,
			 'l': self.vertex_list_l,
			 'm': self.vertex_list_m,
			 'n': self.vertex_list_n,
			 'o': self.vertex_list_o,
			 'p': self.vertex_list_p,
			 'q': self.vertex_list_q,
			 'r': self.vertex_list_r,
			 's': self.vertex_list_s,
			 't': self.vertex_list_t,
			 'u': self.vertex_list_u,
			 'v': self.vertex_list_v,
			 'w': self.vertex_list_w,
			 'x': self.vertex_list_x,
			 'y': self.vertex_list_y,
			 'z': self.vertex_list_z, 
			 
			 '0': self.vertex_list_0,
			 '1': self.vertex_list_1,
			 '2': self.vertex_list_2,
			 '3': self.vertex_list_3,
			 '4': self.vertex_list_4,
			 '5': self.vertex_list_5,
			 '6': self.vertex_list_6,
			 '7': self.vertex_list_7,
			 '8': self.vertex_list_8,
			 '9': self.vertex_list_9,
			 
			 ':': self.vertex_list_colons,
			 ' ': False
			 }
			 
  
  def change_color(self,color) :
    self.color=color
  
  def change_width(self,line_width) :
    self.line_width=line_width
  
  def construct_string(self,string,start_coords) :
    self.string=string
    self.string_to_display=[]
    offset=1.0/self.size
    self.length=0.0
    for v in string :
      char=self.string_look_up.get(v)
      if not char :
	offset += self.size*4
	self.string_to_display.append([])
	self.length += self.size*offset
	continue
      if v == ":" :
	tmp=[]
	for y in char :
	  tmp.append(self._translate(self._translate(y,self.size*offset,0,0),start_coords.wx,start_coords.wy,start_coords.wz))
	self.string_to_display.append(tmp)
	offset += self.size*(2.0/self.size)
	self.length += self.size*offset
	continue
      
      tmp=[]
      for y in char :
	tmp.append(self._translate(self._translate(y,self.size*offset,0,0),start_coords.wx,start_coords.wy,start_coords.wz))
      self.string_to_display.append(tmp)	    
      
      offset += self.size*4
      self.length += self.size*offset
      
        
  def display(self) :
    glEnable(GL_LINE_SMOOTH)
    glHint(GL_LINE_SMOOTH_HINT,GL_NICEST)
    glLineWidth(self.line_width)
    #glEnable(GL_LINE_STIPPLE)
    #glLineStipple( len(self.string_to_display)/8 , len(self.string_to_display)*len(self.string_to_display) ) 
    glLineWidth(self.line_width)
    glColor3ubv(self.color)
    i=0
    for v in self.string_to_display :
      if self.string[i] == ':' :
	glPointSize(int(self.size*6))
	glBegin(GL_POINTS)
	glVertex3fv(v[0].get_vertex())
	glVertex3fv(v[1].get_vertex())
	glEnd()
	i += 1
	continue
 	
      glBegin(GL_LINE_STRIP)
      for y in v :
	glVertex3fv(y.get_vertex())
      glEnd() 
      
      i += 1            
  
  def _translate(self,pt,value_x,value_y,value_z) :
    x=(value_x+pt.wx)
    y=(value_y+pt.wy)
    z=(value_z+pt.wz)
    
    return Vertex(x,y,z)
  
class Start_Anim(Main_Cube_Base_Class) :
  
  def __init__(self,color,number_of_childs=3,max_childs=5,total_childs=25,size=3) :
     
    self.self_color=color 
    
    self.size=size
    
    self.number_of_childs_to_add=number_of_childs
    self.max_childs=max_childs
    self.total_childs=total_childs
    self.death_child_counter=0
    
    self.cubes=[]
    
    self.load_shared_datas() 
    
    
  def load_shared_datas(self) :
    self.initialize_shared_datas()  # Method inherit from the base class. 
    ############################################################################################################################################################################
    # Sequences of Quads:
    #
    # Front-Back: 
    self.quad_middles_front_back=[self.middle_front_back_at_left_down,self.middle_front_back_at_right_down,self.middle_front_back_at_right_up,self.middle_front_back_at_left_up] 
    # Left-Right:
    self.quad_middles_left_right=[self.middle_left_right_at_front_down,self.middle_left_right_at_back_down,self.middle_left_right_at_back_up,self.middle_left_right_at_front_up] 
    # Up-Down: 
    self.quad_middles_up_down=[self.middle_up_down_at_front_right,self.middle_up_down_at_front_left,self.middle_up_down_at_back_left,self.middle_up_down_at_back_right] 
    ############################################################################################################################################################################
    
    
    ############################################################################################################################################################################
    # Sequences of Lines:
    #
    # Right-Left centers:
    self.centers_left_right = [self.center_left_side,self.center_right_side,self.center_left_side,self.center_right_side]
    # Front-Back:
    self.centers_front_back = [self.center_front_side,self.center_back_side,self.center_front_side,self.center_back_side]
    # Up-Down:
    self.centers_up_down    = [self.center_down_side,self.center_up_side,self.center_down_side,self.center_up_side]
    ############################################################################################################################################################################
    
    ######################################################################################################################################################################################
    # Drawing containers:                                                                                                                                                             
    #                                                                                                                                                                                 
    self.quads        = [(self.front,self.back),(self.center_left_side,self.center_right_side),(self.center_front_side,self.center_back_side),(self.center_down_side,self.center_up_side)]
    self.cross_quads  = [self.quad_middles_front_back,self.quad_middles_left_right,self.quad_middles_up_down]                                                                         
    ######################################################################################################################################################################################

  def init_inside_room(self,) :
    ''' Childs cubes generating at start. '''
    
    self.add_ia_cube((randint(127,255),randint(127,255),randint(127,255)),randint(1,6)) 
    
    return
  
  
  def add_ia_cube(self,color,life_points) :
    ''' Function to random the start position and trajectory from an added child cube. '''
    
    # Instantiate an child cube with random settings.
    cube_ia=Cube_IA(color,life_points,False) 
       
    # Append the new child cube to the main cubes childs list.
    self.cubes.append(cube_ia)
    
  def manage_inside_room(self) :
    ''' Main child cubes moving function. '''
    for v in self.cubes :
      v.slide()  
    
  
  def display(self) :
    ''' Main cube display function. '''
    
    glColor3ubv(self.self_color)

    set_lined(self.size)
      
      
    for v in self.quads[0] :
      # Display 2 big quads.
      glBegin(GL_POLYGON)
      for y in v :
	glVertex3fv(y.get_vertex())
      glEnd()  
    
    
    for y in range(0,4) :
      # Join the 2 big quads.
      glBegin(GL_LINES)
      glVertex3fv(self.quads[0][0][y].get_vertex())
      glVertex3fv(self.quads[0][1][y].get_vertex())
      glEnd()
    
    
    
    for v in self.cross_quads :
      # Display the cross from every cube side.
      glBegin(GL_QUADS)
      for y in v :
	glVertex3fv(y.get_vertex())
      glEnd()
    
    i=-1
    while i > -4 :
      # Join the the centers from every side.
      glBegin(GL_LINES)
      for v in self.quads[i] :
	glVertex3fv(v.get_vertex())
      glEnd()
      i -= 1  

class Player_thumbnail() :
  def __init__(self,size,color,line_size=3) :
    
    # Generate base player dislaying datas: triangles, 
    #                                       player middle point computing points list,
    #                                       player cube circonscribe points list for colliding with an cube.
    self.player_display,self.player_center_points,self.collide_center_points=gen_player(1.75)
    #########################################################################################################
    
    
    
    self.self_color=color         # Displaying color data.
    self.size=line_size           # Displaying size data. 
    
    # Set the player in the right start direction.
    m=Matrix()
    m.rotate_y(90)
    self._update_pos(m)
    
    # Translate the player thumbnail to his final place. 
    m=Matrix()
    m.translate((-25.5,-19.5,0.0))
    self._update_pos(m)
    
  def turn_around(self) :
    ''' Function to rotate the player thumbnail '''
    
    # Compute the player thumbnail center.
    self.middle_point=compute_middle_point(self.player_center_points)
    
    # Translate the player thumbnail center to the origin and perform 
    # an rotation on his own Y axe.
    m=Matrix()
    m.translate((-self.middle_point.wx,-self.middle_point.wy,-self.middle_point.wz))
    m.rotate_y(1)
    m.translate((self.middle_point.wx,self.middle_point.wy,-self.middle_point.wz))
    
    self._update_pos(m)   
  
  def display(self) :
    ''' Player display function. '''
    
    glColor3ubv(self.self_color)
    set_lined(self.size)
    
    for v in self.player_display :
      glBegin(GL_TRIANGLES)
      for y in v :
	glVertex3fv(y.get_vertex())
      glEnd()
  
  def _update_pos(self,matrix1) :
    ''' Player slide and colliding positions updates. '''
    
    tmp=[]
    for v in self.player_center_points :
      tmp.append(matrix1.mult_vertex(v))
      
    self.player_center_points=tmp
    
    self._update_coords()
	
      
    
    
    
    
  
  def _update_coords(self) :
    ''' Player display datas updating. '''
    
    pointed=self.player_center_points.pop(0)
    
    self.player_display=[]
    i=0
    while i < len(self.player_center_points) :
      # Loop to build triangles for player pentahedron displaying.
      tmp=[]
      if i+1 != len(self.player_center_points) : 
        tmp.append(self.player_center_points[i])
        tmp.append(self.player_center_points[i+1])
        tmp.append(pointed)
      else :
        tmp.append(self.player_center_points[i])
        tmp.append(self.player_center_points[0])
        tmp.append(pointed)  
      self.player_display.append(tmp)   
      i += 1 
    
    self.player_center_points.insert(0,pointed)
    
class Cube_thumbnail() :
  
  def __init__(self,color,mode,pos_matrix) :
    
    side_length=1.75 # Half player thumbnail size.
    
    self.self_vertex_list=[# Front face vertices:
                           Vertex(-side_length/2.,-side_length/2., -side_length/2.), Vertex(side_length/2.,-side_length/2., -side_length/2.), Vertex(side_length/2., side_length/2., -side_length/2.), Vertex(-side_length/2., side_length/2., -side_length/2.),
                           # Back face vertices.
                           Vertex(-side_length/2.,-side_length/2.,  side_length/2.), Vertex(side_length/2.,-side_length/2.,  side_length/2.), Vertex(side_length/2., side_length/2.,  side_length/2.), Vertex(-side_length/2., side_length/2.,  side_length/2.)]
    
    
    
    
    
    
    self.size=3  # Line width.
    
    self.self_color=color # Player thumbnail color.
    
   
    self._update_pos(pos_matrix) 
     
    self._update_coords()             # Update the datas needed for cube displaying
    
      
    self.middle_point=compute_middle_point(self.self_vertex_list)
    
    
    self.mode=mode
    
    if self.mode == "filled" :
      self.self_border_color=(127,127,127)
    
  def turn_around(self) :
    ''' Cube thumbnail rotating methode. '''
    
    # Compute the cube thumbnail center.
    self.middle_point=compute_middle_point(self.self_vertex_list)
    
    # Translate the cube thumbnail center to the origin and perform 
    # an rotation on his own Y axe.
    m=Matrix()
    m.translate((-self.middle_point.wx,-self.middle_point.wy,-self.middle_point.wz))
    m.rotate_y(1)
    m.translate((self.middle_point.wx,self.middle_point.wy,-self.middle_point.wz))
    self._update_pos(m)  
  

  def _update_pos(self,matrix) :
    ''' updating the displaying datas. '''
    tmp_1=[]
    
    for v in self.self_vertex_list :
      tmp_1.append(matrix.mult_vertex(v))
    
    
    self.self_vertex_list=tmp_1
    self.middle_point=compute_middle_point(tmp_1)
    
    self._update_coords()
    
  def _update_coords(self) :
    ''' Displaying datas updating function. ''' 
    self.front = [self.self_vertex_list[0],self.self_vertex_list[1],self.self_vertex_list[2],self.self_vertex_list[3]]
    self.back  = [self.self_vertex_list[4],self.self_vertex_list[5],self.self_vertex_list[6],self.self_vertex_list[7]]
    
    
    self.left  = [self.self_vertex_list[7],self.self_vertex_list[4],self.self_vertex_list[0],self.self_vertex_list[3]]
    self.right = [self.self_vertex_list[5],self.self_vertex_list[6],self.self_vertex_list[1],self.self_vertex_list[2]]
    
    self.up    = [self.self_vertex_list[4],self.self_vertex_list[5],self.self_vertex_list[1],self.self_vertex_list[0]]
    self.down  = [self.self_vertex_list[7],self.self_vertex_list[6],self.self_vertex_list[2],self.self_vertex_list[3]]
      
    
    self.quads_lined_base   = [self.front,self.back]
    self.quads_filled_base  = [self.front,self.back,self.left,self.right,self.up,self.down]
    
  def display(self) :
    
    glColor3ubv(self.self_color)
    
    if self.mode == 'lined' :
      set_lined(self.size)
      
      for v in self.quads_lined_base :
	# We display 2 quads.
        glBegin(GL_LINE_LOOP)
        for y in v :
          glVertex3fv(y.get_vertex())
        glEnd()  
      
      
      for y in range(0,4) :
	# We join the 2 quads.
        glBegin(GL_LINES)
        glVertex3fv(self.quads_lined_base[0][y].get_vertex())
        glVertex3fv(self.quads_lined_base[1][y].get_vertex())
        glEnd()
      
            
          
          
    if self.mode == "filled" :
      set_polygon_mode(mode=GL_FILL)
      
      for v in self.quads_filled_base :
	# Displaying every cube face.
        glBegin(GL_QUADS) 
        for y in v :
          glVertex3fv(y.get_vertex())
        glEnd() 
      
      set_lined(self.size)
      glColor3ubv(self.self_border_color)
      for v in self.quads_lined_base :
	# Displaying every cube border.
        glBegin(GL_LINE_LOOP)
        for y in v :
          glVertex3fv(y.get_vertex())
        glEnd()  
      
      
      for y in range(0,4) :
        glBegin(GL_LINES)
        glVertex3fv(self.quads_lined_base[0][y].get_vertex())
        glVertex3fv(self.quads_lined_base[1][y].get_vertex())
        glEnd()      
    
  
def move_player(event) :
  ''' Player control function throught keyboard events given as argument. '''
  
  if event.type == QUIT :
    exit()
      
  elif event.type == KEYDOWN :
    
    if event.key == K_UP :
      # Call the player orientation function.
      player.change_player_direction('U')  
    
    elif event.key == K_DOWN :
      # Call the player orientation function.
      player.change_player_direction('D') 
      
    elif event.key == K_RIGHT :
      # Call the player orientation function.
      player.change_player_direction('R')  
    
    elif event.key == K_LEFT :
      # Call the player orientation function.
      player.change_player_direction('L')    
    
    elif event.key == K_z or event.key == K_w :
      # Main player move callback.
      player.move()
    
    elif event.key == K_s :
      # Main player shooting callback.
      player.shoot()
    
    elif event.key == K_d  and player.is_dead : 
      # Let appear the player in the main cube center.
      player.new_life()
      
    elif event.key == K_SPACE :
      change_game_mode(False)
  

def move_playground(event) :
  
  global c_right,c_up,c_sight # This are the X, Y, Z Localview axe.
  
  
  
  if event.type == QUIT :
    exit()
      
  elif event.type == KEYDOWN :
    if event.key == K_UP :
      # Rotate the playground around the X axe.
      c_up += 0.75           
      camera()  
    
    elif event.key == K_DOWN :
      # Rotate the playground around the X axe.
      c_up -= 0.75
      camera()
      
    elif event.key == K_RIGHT :
      # Rotate the playground around the Z axe.
      c_right += 0.75
      camera()   
    
    elif event.key == K_LEFT :
      # Rotate the playground around the Z axe.
      c_right -= 0.75
      camera()    
    
    elif event.key == K_d :
      # Rotate the playground around the Y axe.
      c_sight += 0.75                              
      camera() 
    
    elif event.key == K_q or event.key == K_a :
      # Rotate the playground around the Y axe.
      c_sight -= 0.75
      camera()
    
    elif event.key == K_z or event.key == K_w :
      # Reset the playground to the start view.
      reset_camera()
    
    elif event.key == K_KP_PLUS :
      # Increment the lines width.
      
      if main_cube.size < 7 :
        main_cube.size += 1
        
        for v in main_cube.cubes :
          v.size += 1
          
        player.size += 1
        
    elif event.key == K_KP_MINUS :
      # Decrement the lines width.
      
      if main_cube.size > 2 :
        main_cube.size -= 1
        
        for v in main_cube.cubes :
          v.size -= 1
          
        player.size -= 1 
        
    elif event.key == K_SPACE :
      change_game_mode(True)
    
    
  

def change_game_mode(mode) :
  ''' Change the display actions mode, even in:
      -) The user play the game.
      -) The user can rotate the main cube and scale the line width.
  ''' 
  
  global game_run
  
  game_run=mode  # Boolean variable for set the mode.
  
  if mode == False :
    # Settings for the display configuration mode.
    sleeper=0.001
    pygame.key.set_repeat(100,10)
    
  elif mode == True :
    # Settings for the game playing mode.
    sleeper=0.15
    pygame.key.set_repeat(False,False)

def camera() :
  ''' Camera settings function. '''
  camera_view=Localview()                     # instantiate an localview
  m=Matrix()                                  # Build the camera position 
  m.translate((0.0,0.0,0.0))                  # Build the camera position 
  m.rotate_vector(c_right,camera_view.right)  # Build the camera position 
  m.rotate_vector(c_up,camera_view.up)        # Build the camera position 
  m.rotate_vector(c_sight,camera_view.sight)  # Build the camera position 
  m.load_hardware()                           # Load matrix as MODELVIEW matrix
  
  
  
def reset_camera() :
  global c_right,c_up,c_sight
  c_right,c_up,c_sight=22.5*1.5, 22.5*1.5,-22.5*1.5  
  camera()
  
class Mainloop() :
  def __init__(self) :
    ''' Settings by game start. '''
    self.at_start=True
    
    # Default settings:
    self.player_lifes=5
    self.number_of_cubes_in_playground=5
    self.total_number_of_cubes=25  
  
  def init_objects(self) :
    ''' Initialise an game with the choosen settings.'''
    
    global main_cube, player, cube, player, game_run, sleeper, c_right, c_up, c_sight
    
    c_right,c_up,c_sight=22.5*1.5, 22.5*1.5,-22.5*1.5
    
    game_run=True
    
    player=Player(1.75,(255,0,0),self.player_lifes)
    
    main_cube=Main_Cube((127,127,127),self.number_of_cubes_in_playground,self.total_number_of_cubes)
    main_cube.init_inside_room()
    
    sleeper=0.15
    
    
    
  

  def change_selection(self,selection,prev_selection) :
    ''' Change settings calllback by game configuration. '''
    
    if selection == 3 :
      # The user is by the start game choice. 
      
      # We highlight the selection line.  
      self.start_game_text.change_color((127,127,127))
      self.start_game_text.change_width(6)
      self.total_number_of_cubes_text.change_color((255,0,0))
      self.total_number_of_cubes_text.change_width(2)
      
    
    elif selection == 2 :
      # The user is by the total number of cubes to shoot down setting line.
      
      # We highlight the line.
      self.total_number_of_cubes_text.change_color((127,127,127))
      self.total_number_of_cubes_text.change_width(3)
      if prev_selection == 3 :
	# We reconfigure the precedent configuration line.
	self.start_game_text.change_color((255,0,0))
	self.start_game_text.change_width(2)
      elif prev_selection == 1 :
	# We reconfigure the precedent configuration line.
	self.cubes_in_playground_text.change_color((255,0,0))
	self.cubes_in_playground_text.change_width(2)
	
    elif selection == 1 :
      # The user is by the number of cubes simultaneaous in the playground setting line.
      
      # We highlight the line.
      self.cubes_in_playground_text.change_color((127,127,127))
      self.cubes_in_playground_text.change_width(3)
      if prev_selection == 2 :
	# We reconfigure the precedent configuration line.
	self.total_number_of_cubes_text.change_color((255,0,0))
	self.total_number_of_cubes_text.change_width(2)
      elif prev_selection == 0 :
	# We reconfigure the precedent configuration line.
	self.number_of_life_text.change_color((255,0,0))
	self.number_of_life_text.change_width(2)
	
    elif selection == 0 :
      # The user is by the number of player lifes setting line.
      
      # We highlight the line.
      self.number_of_life_text.change_color((127,127,127))
      self.number_of_life_text.change_width(3)
      if prev_selection == 1 :
	# We reconfigure the precedent configuration line.
	self.cubes_in_playground_text.change_color((255,0,0))
	self.cubes_in_playground_text.change_width(2)
      elif prev_selection == 3 :
	# We reconfigure the precedent configuration line.
	self.start_game_text.change_color((255,0,0))
	self.cubes_in_playground_text.change_width(5)
	
      
      
  
  def configure_game(self) :
    ''' Method to let the user configurate his game with the wanted settings. '''
    
    self.at_start=False
    
    main_matrix.load_hardware()
    
    # Creating decorative player and cubes thumbnails.
    player=Player_thumbnail(1.75,(255,0,0))
    
    m=Matrix()
    m.translate((-25.5,-9.25,0.0))
    cube1= Cube_thumbnail((255,0,0),"lined",m)
    
    m=Matrix()
    m.translate((-25.5,0.75,0.0))
    cube2= Cube_thumbnail((255,0,0),"filled",m)
    
    
    setting_selected=3 # At begin the choice is to start the game.
    
    
    
    # Creating the number of life setting configuration text.
    self.number_of_life_text=Text(size=0.65,line_width=2,color=(255,0,0))
    self.number_of_life_text.construct_string("player number of lifes  : %s" % str(self.player_lifes).zfill(2),Vertex(-22.5,-20.0,0.0))
    
    # Creating the maximal number of cubes in the playground setting configuration text.
    self.cubes_in_playground_text=Text(size=0.65,line_width=2,color=(255,0,0))
    self.cubes_in_playground_text.construct_string("number of cubes simultan: %s" % str(self.number_of_cubes_in_playground).zfill(2) ,Vertex(-22.5,-10.0,0.0))
    
    # Creating the total number of cubes tro shoot down setting configuration text.
    self.total_number_of_cubes_text=Text(size=0.65,line_width=2,color=(255,0,0))
    self.total_number_of_cubes_text.construct_string("total number of cubes   : %s" % str(self.total_number_of_cubes).zfill(2),Vertex(-22.5,0.0,0.0))
    
    # Creating the game starting configuration text.
    self.start_game_text=Text(size=0.95,line_width=5,color=(127,127,127))
    self.start_game_text.construct_string("start game",Vertex(-17.5,10.0,0.0))
    
    while True :
      
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT) 
	
      
      player.turn_around()
      player.display()
      
      cube1.turn_around()
      cube1.display()
      
      cube2.turn_around()
      cube2.display()
      
      self.number_of_life_text.display()
      
      self.cubes_in_playground_text.display()
      
      self.total_number_of_cubes_text.display()
      
      self.start_game_text.display()
      
      for event in pygame.event.get() :
	
	if event.type == QUIT :
	  exit()
	    
	if event.type == KEYDOWN :
	  
	  if event.key == K_UP :
	    if setting_selected > 0 and setting_selected <= 3 :
	      # Change configuration settings line. 
	      prev_selection=setting_selected
	      setting_selected -= 1
	      self.change_selection(setting_selected,prev_selection)
	   	    
	  elif event.key == K_DOWN :
	    # Change configuration settings line. 
	    if setting_selected >= 0 and setting_selected < 3 :
	      prev_selection=setting_selected
	      setting_selected += 1
	      self.change_selection(setting_selected,prev_selection)
	      
	  if event.key == K_RIGHT and setting_selected == 0 :
	    # Change the value of players life setting by increment it. 
	    if self.player_lifes < 99 and self.player_lifes >= 1 :
	      self.player_lifes += 1
	      self.number_of_life_text=Text(size=0.65,line_width=3,color=(127,127,127))
	      self.number_of_life_text.construct_string("player number of lifes  : %s" % str(self.player_lifes).zfill(2),Vertex(-22.5,-20.0,0.0))
	      
	  elif event.key == K_LEFT and setting_selected == 0 :
	    # Change the value of players life setting by decrement it. 
	    if self.player_lifes <= 99 and self.player_lifes > 1 :
	      self.player_lifes -= 1
	      self.number_of_life_text=Text(size=0.65,line_width=3,color=(127,127,127))
	      self.number_of_life_text.construct_string("player number of lifes  : %s" % str(self.player_lifes).zfill(2),Vertex(-22.5,-20.0,0.0)) 
	      
	  if event.key == K_RIGHT and setting_selected == 1 :
	    # Change the maximal value of cubes in the playground setting by increment it. 
	    if self.number_of_cubes_in_playground < 99 and self.number_of_cubes_in_playground >= 1 :
	      self.number_of_cubes_in_playground += 1
	      self.cubes_in_playground_text=Text(size=0.65,line_width=3,color=(127,127,127))
	      self.cubes_in_playground_text.construct_string("number of cubes simultan: %s" % str(self.number_of_cubes_in_playground).zfill(2),Vertex(-22.5,-10.0,0.0))
	  
	  elif event.key == K_LEFT and setting_selected == 1 :
	    # Change the maximal value of cubes in the playground setting by decrement it. 
	    if self.number_of_cubes_in_playground <= 99 and self.number_of_cubes_in_playground > 1 :
	      self.number_of_cubes_in_playground -= 1
	      self.cubes_in_playground_text=Text(size=0.65,line_width=3,color=(127,127,127))
	      self.cubes_in_playground_text.construct_string("number of cubes simultan: %s" % str(self.number_of_cubes_in_playground).zfill(2),Vertex(-22.5,-10.0,0.0)) 
	      
	  if event.key == K_RIGHT and setting_selected == 2 :
	    # Change the value of cubes total to shoot down by increment it. 
	    if self.total_number_of_cubes < 99 and self.total_number_of_cubes >= 1 :
	      self.total_number_of_cubes += 1
	      self.total_number_of_cubes_text=Text(size=0.65,line_width=3,color=(127,127,127))
	      self.total_number_of_cubes_text.construct_string("total number of cubes   : %s" % str(self.total_number_of_cubes).zfill(2),Vertex(-22.5,0.0,0.0))
	  
	  elif event.key == K_LEFT and setting_selected == 2 :
	    # Change the value of cubes total to shoot down by decrement it. 
	    if self.total_number_of_cubes <= 99 and self.total_number_of_cubes > 1 and self.total_number_of_cubes > self.number_of_cubes_in_playground :
	      self.total_number_of_cubes -= 1
	      self.total_number_of_cubes_text=Text(size=0.65,line_width=3,color=(127,127,127))
	      self.total_number_of_cubes_text.construct_string("total number of cubes   : %s" % str(self.total_number_of_cubes).zfill(2),Vertex(-22.5,0.0,0.0))       
	  
	  elif event.key == K_RETURN and setting_selected == 3 :
	    # Start the game. 
	    return
	  
      pygame.display.flip()    
      sleep(0.015)
      
    
 
  def start_anim(self) :
    ''' Game start animation '''
    global c_right,c_up,c_sight
    c_right,c_up,c_sight=0,0,0
    
    main_cube=Start_Anim((127,127,127))
    
    if self.at_start :
      # The user as start the programme.
      
      # Load and play the background music.
      pygame.mixer.music.load("/usr/share/cube-hunter/Sound/Music/cuber hunter music.mp3")
      pygame.mixer.music.play(-1)
      
      # Start Localview axes value.
      axes_value=22.5*1.5/45.0  
      
      # Before the Main Cube rotation animation no child cubes are displayed.
      cubes_display_init=False
      
      # The titel line is not displayed on the top of the display.
      display_titel=False
      
      # Animation speed at begin.
      sleeper=0.15*1.65
      
      # Main cube displaying position at start.
      camera_view=Localview()                     # instantiate an localview
      m=Matrix()                                  # Build the camera position 
      m.translate((0.0,0.0,0.0))                  # Build the camera position 
      m.rotate_vector(c_right,camera_view.right)  # Build the camera position 
      m.rotate_vector(c_up,camera_view.up)        # Build the camera position 
      m.rotate_vector(c_sight,camera_view.sight)  # Build the camera position 
      m.load_hardware() 
      
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT) 
      
      
	
      camera() 
      
      pygame.display.flip() 
      
      sleep(3)
    
    else :
      # The user has soon sea the start animation.
      
      # Load and play the background music.
      pygame.mixer.music.load("/usr/share/cube-hunter/Sound/Music/cube hunter music 2.wav")
      pygame.mixer.music.play(-1)
      
      # Start Localview axes value.
      axes_value=22.5*1.5
      
      # After the Main Cube rotation animation the child cubes are displayed.
      cubes_display_init=True
      
      # The titel line is displayed on the top of the display.
      display_titel=True
      
      # Animation speed.
      sleeper=0.15
      
      # Main cube displaying position.
      camera_view=Localview()                           # instantiate an localview
      m=Matrix()                                        # Build the camera position 
      m.translate((0.0,0.0,0.0))                        # Build the camera position 
      m.rotate_vector(axes_value,camera_view.right)     # Build the camera position 
      m.rotate_vector(axes_value,camera_view.up)        # Build the camera position 
      m.rotate_vector(axes_value,camera_view.sight)     # Build the camera position 
      m.load_hardware() 
      
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT) 
      
      
	
      camera()
      pygame.display.flip() 
      
    
    player_slide=False
    
    incr=axes_value
    
    anim_main_cube=True
    
    
    add_one_cube=False
    
    wait_add_a_cube_counter=1
    
    # Contruct the Text strings object.
    titel_game_part_1=Text(size=1.0,line_width=7,color=(255,0,0))
    titel_game_part_1.construct_string("cube",Vertex(-19.5,-2.0,0.0))
    
    titel_game_part_2=Text(size=1.0,line_width=7,color=(255,0,0))
    titel_game_part_2.construct_string("hunter",Vertex(-23.5,4.0,0.0))
    
    titel_game=Text(size=1.05,line_width=7,color=(255,0,0))
    titel_game.construct_string("cube hunter",Vertex(-22.0,-28.5,0.0))
    
    press_enter_to_config=Text(size=0.65,line_width=2,color=(255,0,0))
    press_enter_to_config.construct_string("press enter to configure game",Vertex(-25.0,27.5,0.0))
    
    # Press enter to configure game Text blinking time.
    blink_counter=1
    
    while True :
      
      if blink_counter == 11 :
	blink_counter=1
      
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT) 
      
      if display_titel :
	# The titel line is displayed on the top of the display.
	
	glPushMatrix() # Push a matrix to display the Text independant of the Main_Cube position. 
	
	# Settings the Text diplaying position.
	camera_view=Localview()                     # instantiate an localview
	m=Matrix()                                  # Build the camera position 
	m.translate((0.0,0.0,0.0))                  # Build the camera position 
	m.rotate_vector(c_right,camera_view.right)  # Build the camera position 
	m.rotate_vector(c_up,camera_view.up)        # Build the camera position 
	m.rotate_vector(c_sight,camera_view.sight)  # Build the camera position 
	m.load_hardware() 
	
	titel_game.display()
	
	press_enter_to_config.display()
	
	glPopMatrix()
	
	blink_counter += 1
      
      main_cube.manage_inside_room()	  
      main_cube.display()
      
      
      
      
	
      if axes_value < 22.5*1.5 and anim_main_cube :
	# The Main_Cube rotating animation is in action.
	
	camera_view=Localview()                        # instantiate an localview
	m=Matrix()                                     # Build the camera position 
	m.translate((0.0,0.0,0.0))                     # Build the camera position 
	m.rotate_vector(axes_value,camera_view.right)  # Build the camera position 
	m.rotate_vector(axes_value,camera_view.up)     # Build the camera position 
	m.rotate_vector(axes_value,camera_view.sight)  # Build the camera position 
	m.load_hardware() 
	
	titel_game_part_1.display()
	titel_game_part_2.display()
	
	
	  
      elif axes_value == 22.5*1.5 and not add_one_cube :
	# The Main_Cube rotating animation is finish.
	 
	if self.at_start :
	  sleep(3)
	
	# Settings for childs cubes animation.
	cubes_display_init=True
	add_one_cube=True
	anim_main_cube=False
	wait_add_a_cube_counter=50
	
	camera_view=Localview()                           # instantiate an localview
	m=Matrix()                                        # Build the camera position 
	m.translate((0.0,0.0,0.0))                        # Build the camera position 
	m.rotate_vector(axes_value,camera_view.right)     # Build the camera position 
	m.rotate_vector(axes_value,camera_view.up)        # Build the camera position 
	m.rotate_vector(axes_value,camera_view.sight)     # Build the camera position 
	m.load_hardware() 
    
      if cubes_display_init and wait_add_a_cube_counter % 50 == 0 and wait_add_a_cube_counter < 450 :
	main_cube.init_inside_room()
	sleeper=0.15
	cubes_display_init=False
	display_titel=True
      
      
      
      for event in pygame.event.get() :
	if event.type == QUIT :
	  exit()
	  
	if event.type == KEYDOWN :
	  if event.key == K_ESCAPE :
	    exit()
	  elif event.key == K_RETURN :
	    return
	    
	    
      if add_one_cube :
	wait_add_a_cube_counter += 1
	cubes_display_init=True
      
      if wait_add_a_cube_counter > 450 :
	wait_add_a_cube_counter=451
      
      if not axes_value == 22.5*1.5 + incr :
	axes_value += incr
      
      
	
      if blink_counter == 5 :
	press_enter_to_config.change_color((127,127,127))
      elif blink_counter == 10 :
	press_enter_to_config.change_color((255,0,0))  
      
      pygame.display.flip()    
      sleep(sleeper)
  
 
  def game_mainloop(self) :
    
    global game_run
    
    pygame.mixer.music.stop()
    
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT) 
    
    
    camera()
    pygame.display.flip() 
    
    player_slide=False
    
    
    while True :
      
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT) 
      
      
      if game_run :
      
	if player.is_shooting :
	  
	  for v in player.shoot_balls :
	    v.display()
	    v.display()
	    v.display()
	    v.display()
	    
	main_cube.display()
	main_cube.manage_inside_room()
	
	for v in range(0,2) :
	  player_slide=True
	  player.slide()
	  player.collide_check()
	  
	player.display()
	
	if not player_slide :
	  player.collide_check()
	else :
	  player_slide=False
	
	
      else :
	main_cube.display()
	for v in main_cube.cubes :
	  v.display()
	  
	player.display()
      
      for event in pygame.event.get() :
	if event.type == USEREVENT :
	  if event.code == "end_game" :
	    return
	  
	if not game_run :
	  move_playground(event)
	else :
	  move_player(event)
      
      pygame.display.flip()    
      sleep(sleeper)
  

    
    
    
 
  def run(self) :
    ''' The user launch the programme or has finish a game. '''
    
    self.start_anim()      # We display the start animation.
    
    self.configure_game()  # The user can configure his game settings.
    
    self.init_objects()    # Initialise the game with the user settings.
    
    self.game_mainloop()   # Game mainloop
    
    self.run()             # Game finish we recall this method but the configuration inhibit the Main_Cube rotation animation.
  
 
def init_GL() :
  
  glEnable(GL_DEPTH_TEST)          # Enable z-buffer  
  
  
  glClearColor(0.0, 0.0, 0.0, 0.0) # Define clear color [0.0-1.0]
  
  glShadeModel(GL_FLAT)            # Define lines as polygon instead of full polygon: GL_SMOOTH
  
  
  
def resize(width,height) :
  
  fov_angle     = 60.0
  distance_near = 2.0
  distance_far  = 1000.0 
  
  glViewport(0, 0, int(width), int(height))
  glMatrixMode(GL_PROJECTION)
  
  glLoadIdentity()
  gluPerspective( fov_angle, float(width)/float(height),distance_near, distance_far )
  
  glMatrixMode(GL_PROJECTION)
  glLoadIdentity()
  glOrtho(-30.0,30.0,30.,-30.0, -30.0,30.0)
  glMatrixMode(GL_MODELVIEW)
  glLoadIdentity()
   
 
def main():
  global screen,width,height,Tk,main_matrix
  
  screen_config=Tk()
  
  height = screen_config.winfo_screenheight()
  width  = height * (800./600.)
  
  del(screen_config)
  del(Tk)
  
  sizer=0.80
  
  
  width  = width  * sizer
  height = height * sizer 
  
  pygame.init()
  screen=pygame.display.set_mode((int(width),int(height)), HWSURFACE | OPENGL | DOUBLEBUF,24)
  pygame.display.set_caption("cube hunter")
  
  glutInit()
  
  resize(width,height)
  init_GL()
  
  
  main_matrix=Matrix()
  
  main=Mainloop()
  main.run()
  
if __name__ == "__main__" :
  main()  