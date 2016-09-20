#ifndef FLATPRIM_H
#define FLATPRIM_H

#include "SceneData.h"
#include <vector>

struct FlatPrim
{
    std::vector<ScenePrimitive*> prims;
    Matrix matrix;
    Matrix inverted;
};

class PrimList
{
public:
    PrimList(SceneNode *rootNode) {
        //DFS on rootNode to populate vector of children
        Matrix identity;
        populate(rootNode, identity);
    };

    ~PrimList() {
        //Loop through children and delete each
        for(std::vector<FlatPrim*>::iterator it = prims.begin(); it != prims.end(); ++it) {
            delete (*it);
        }
    };

    std::vector<FlatPrim*> prims;
private:
    void populate(SceneNode *root, Matrix prev) {
        if (root == NULL) {
            return;
        }
        Matrix current = getTransformMat(root->transformations, prev);
        Matrix inverted = invert(current);
        if (!root->primitives.empty()) {
            FlatPrim *node = new FlatPrim;
            node->matrix = current;
            node->inverted = inverted;
            node->prims = root->primitives;
            prims.push_back(node);
        }
        for(std::vector<SceneNode*>::iterator it = root->children.begin(); it != root->children.end(); ++it) {
            populate(*it, current);
        }
    };

    Matrix getTransformMat(std::vector<SceneTransformation*> &children, Matrix prev) {
        Matrix rVal = prev;
        for(std::vector<SceneTransformation*>::iterator it = children.begin(); it != children.end(); ++it) {
            Matrix m;
            switch((*it)->type) {
                case TRANSFORMATION_TRANSLATE:
                    m = trans_mat((*it)->translate);
                    break;
                case TRANSFORMATION_SCALE:
                    m = scale_mat((*it)->scale);
                    break;
                case TRANSFORMATION_ROTATE:
                    m = rot_mat((*it)->rotate, (*it)->angle);
                    break;
                case TRANSFORMATION_MATRIX:
                    m = (*it)->matrix;
                    break;
                default:
                    std::exit(0);
            }
            rVal = rVal * m;
        }
        return rVal;
    }
};

#endif
