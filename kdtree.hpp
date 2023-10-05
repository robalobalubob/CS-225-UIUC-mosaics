/**
 * @file kdtree.cpp
 * Implementation of KDTree class.
 */

#include <utility>
#include <algorithm>

using namespace std;

template <int Dim>
bool KDTree<Dim>::smallerDimVal(const Point<Dim>& first,
                                const Point<Dim>& second, int curDim) const
{
    /**
     * @todo Implement this function!
     */
    if (curDim < 0 || curDim >= Dim) {
      return false;
    }
    if (first[curDim] == second[curDim]) {
      return first < second;
    }
    return (first[curDim] < second[curDim]);
}

template <int Dim>
bool KDTree<Dim>::shouldReplace(const Point<Dim>& target,
                                const Point<Dim>& currentBest,
                                const Point<Dim>& potential) const
{
    /**
     * @todo Implement this function!
     */
    float dist = 0;
    float potDist = 0;
    for (unsigned i = 0; i < Dim; i++) {
      dist += ((target[i]-currentBest[i]) * (target[i]-currentBest[i]));
      potDist += ((target[i]-potential[i])*(target[i]-potential[i]));
    }
    if (dist == potDist) {
      return potential < currentBest;
    }
    return potDist < dist;
}


template<int Dim>
Point<Dim>& KDTree<Dim>::qSelect(vector<Point<Dim>>& list, int dimen, unsigned left, unsigned right, unsigned i) {
  if (left == right) {
    return list[left];
  }
  unsigned piv = i;
  
  Point<Dim> piv_val = list[piv];
  Point<Dim> t = list[piv];
  list[piv] = list[right];
  list[right] = t;
  unsigned ind = left;
  for(unsigned j = left; j < right; j++) {
    if(smallerDimVal(list[j], piv_val, dimen)) {
      t = list[ind];
      list[ind] = list[j];
      list[j] = t;
      ind++;
    }
  }
  t = list[ind];
  list[ind] = list[right];
  list[right] = t;
  piv = ind;

  if (i == piv) {
    return list[piv];
  } else if (i < piv) {
    return qSelect(list, dimen, left, piv - 1, i);
  } else {
    return qSelect(list, dimen, piv + 1, right, i);
  }

}

template<int Dim>
typename KDTree<Dim>::KDTreeNode* KDTree<Dim>::constHelp(vector<Point<Dim>>& points, int dimen, unsigned left, unsigned right) {
  if(points.empty() || left < 0 || right>=points.size() || right < 0 || left >= points.size()) {
    return NULL;
  }
  if (left > right) {
    return NULL;
  }
  unsigned med_ind = (left + right) / 2;
  KDTreeNode* subroot = new KDTreeNode(qSelect(points, dimen%Dim, left, right, med_ind));
  size += 1;
  dimen++;
  subroot->left = constHelp(points, dimen, left, med_ind - 1);
  subroot->right = constHelp(points, dimen, med_ind + 1, right);
  return subroot;
}

template <int Dim>
KDTree<Dim>::KDTree(const vector<Point<Dim>>& newPoints)
{
    /**
     * @todo Implement this function!
     */
    size = 0;
    vector<Point<Dim>> points;
    points.assign(newPoints.begin(), newPoints.end());
    root = constHelp(points, 0, 0, points.size() - 1);

}

template <int Dim>
KDTree<Dim>::KDTree(const KDTree<Dim>& other) {
  /**
   * @todo Implement this function!
   */
  copy(this->root, other->root);
  size = other.size;
}

template <int Dim>
const KDTree<Dim>& KDTree<Dim>::operator=(const KDTree<Dim>& rhs) {
  /**
   * @todo Implement this function!
   */
  if (this != NULL) {
    clear(root);
  }
  copy(this->root, rhs->root);
  size = rhs.size;
  return *this;
}

template <int Dim>
KDTree<Dim>::~KDTree() {
  /**
   * @todo Implement this function!
   */
  clear(root);
}

template <int Dim>
Point<Dim> KDTree<Dim>::findNearestNeighbor(const Point<Dim>& query) const
{
    /**
     * @todo Implement this function!
     */

    return findNearestNeighbor(root, query, 0);
}

template <int Dim>
Point<Dim> KDTree<Dim>::findNearestNeighbor(KDTreeNode* subroot, const Point<Dim>& query, int dim) const {
  Point<Dim> curBest = subroot->point;
  bool yoyo;
  if(subroot->left == NULL && subroot->right == NULL) {
    return subroot->point;
  }
  if(smallerDimVal(query,subroot->point, dim)) {
    if (subroot->left == NULL) {
      curBest = findNearestNeighbor(subroot->right, query, (dim+1) % Dim);
    
    } else {
      curBest = findNearestNeighbor(subroot->left, query, (dim+1) % Dim);
    }
    yoyo = true;
  } else {
    if (subroot->right == NULL) {
      curBest = findNearestNeighbor(subroot->left, query, (dim+1) % Dim);
    } else {
      curBest = findNearestNeighbor(subroot->right, query, (dim+1) % Dim);
    }
    yoyo = false;
  }
  if (shouldReplace(query, curBest, subroot->point)) {
    curBest = subroot->point;
  }
  float rad = 0;
  for (int i = 0; i < Dim; i++) {
    rad += (query[i] - curBest[i]) * (query[i] - curBest[i]);
  }
  float dist = subroot->point[dim] - query[dim];
  dist = dist * dist;
  if (dist <= rad) {

    KDTreeNode* check;
    if (yoyo) {
      check = subroot->right;
    } else {
      check = subroot->left;
    }
    if (check != NULL) {
      Point<Dim> bester = findNearestNeighbor(check, query, (dim+1)%Dim);
      if (shouldReplace(query, curBest, bester)) {
        curBest = bester;
      }
    }
  }
  return curBest;
}

template <int Dim>
void KDTree<Dim>::clear(KDTreeNode* subroot) {
  if (subroot == NULL) {
    return;
  }
  if (subroot->left != NULL) {
    clear(subroot->left);
  }
  if(subroot->right != NULL) {
    clear(subroot->right);
  }
  delete subroot;
  subroot = NULL;
}

template <int Dim>
void KDTree<Dim>::copy(KDTreeNode* subroot, KDTreeNode* root2) {
  subroot = new KDTreeNode();
  subroot->point = root2->point;
  copy(subroot->left, root2->left);
  copy(subroot->right, root2->right);
}