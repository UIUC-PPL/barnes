/*
 * ForceData:
 * Has fields to store the moments of a tree node. This includes the 
 * center of mass, total mass and radius of node.
 */
class ForceData {
  OrientedBox<double> box = new OrientedBox<double>();
  MultipoleMoments moments = new MultipoleMoments();

  OrientedBox<double> getBox(){
    return box;
  }

  MultipoleMoments getMoments(){
    return moments;
  }
};
