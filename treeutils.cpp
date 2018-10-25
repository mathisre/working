//#include <stdio.h>
//#include <math.h>
#include <string>

#include "treeutils.h"
using namespace std;

devtree *devtree::hasdev(devtree *root, int val1)
{
  if (root == NULL) return NULL;
  //printf("hasdev() %d %p %p %p\n",val1,root,root->less, root->more);
  if (val1<root->val)
    {
      if (root->less==NULL) return NULL;
      else return hasdev(root->less,val1);
    }
  else if (val1>root->val)
    {
      if (root->more==NULL) return NULL;
      else return hasdev(root->more,val1);
    }
  else //val1==root->val
    {
      return root;
    }
}

devtree* devtree::remove(devtree *dt)
{
  if (dt == NULL) return NULL;

  devtree *minofmore, *minofmoreowner, *retval, *parent;

  //printf("remove %d\n",dt->val);

  parent = dt->parent;

  if (dt->more==NULL)
    {
      //printf("more==NULL\n");
      retval=dt->less;
      if (retval!=NULL) retval->parent=parent;
    }
  else
    {
      //printf("more!=NULL\n");
      minofmore = findmin(dt->more);
      minofmoreowner = minofmore->parent;
      retval = minofmore;
      retval->parent=dt->parent;
      retval->less=dt->less;
      if (retval->less!=NULL) retval->less->parent=retval;

      if (minofmoreowner != dt)
	{
	  //printf("minofmoreowner != dt\n");
	  dt->more->parent=retval;
	  minofmoreowner->less=retval->more;
	  if (retval->more!=NULL) retval->more->parent=minofmoreowner;
	  retval->more = dt->more;
	}
    }
  if (parent!=NULL)
    {
      //printf("parent != NULL\n");
      if (parent->less==dt) parent->less=retval;
      else if (parent->more==dt) parent->more=retval;
    }
  delete dt;
  return retval;
}

devtree *devtree::add(devtree *dt,int val1)
{
  if (dt == NULL) return new devtree(val1);
  devtree *dv2;

  if (val1>=dt->val)
    {
      if (dt->more==NULL)
	{
	  dv2=new devtree(val1);
	  dt->more=dv2;
	  dv2->parent=dt;
	  return dv2;
	}
      else return add(dt->more,val1);
    }
  else
    {
      if (dt->less==NULL)
	{
	  dv2=new devtree(val1);
	  dt->less=dv2;
	  dv2->parent=dt;
	  return dv2;
	}
      else return add(dt->less,val1);
    }
}

devtree* devtree::findmin(devtree *root)
{
  if (root->less==NULL) return root;
  else return findmin(root->less);
}

devtree *devtree::successor(devtree *start)
{
  if (start == NULL) return NULL;
  devtree *parent;
  if (start->more!=NULL)
    return findmin(start->more);
  else
    {
      parent = start->parent;
      while (parent!=NULL && parent->more==start)
	{
	  start=parent;
	  parent=start->parent;
	}
      return parent;
    }
}

bool devtree::equal(devtree *r1, devtree *r2)
{
  r1=findmin(r1);
  r2=findmin(r2);
  while ((r1!=NULL) && (r2!=NULL))
    {
      if (r1->val!=r2->val) return false;
      r1 = successor(r1);
      r2 = successor(r2);
    }
  if ((r1==NULL) && (r2==NULL)) return true;
  return false;
}

int *devtree::vals(devtree *root, int size)
{
  int i;
  int * retval;

  if (size ==0) return NULL;
  retval = new int[size];
  i=0;
  root = findmin(root);
  while ((root!=NULL) && (i<size))
    {
      retval[i]=root->val;
      root = successor(root);
      i++;
    }
  return retval;
}

void devtree::deltree(devtree *root)
{
  if (root==NULL) return;
  //  printf("deleting val: %d %p %p %p\n",root->val,root,root->less,root->more);
  if (root->less!=NULL) deltree(root->less);
  if (root->more!=NULL) deltree(root->more);
  delete root;
}
