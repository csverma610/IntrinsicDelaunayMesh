#include "IntrinsicDelaunay.h"

using namespace std;
#include <fstream>

////////////////////////////////////////////////////////////////////////////////

void IntrinsicDelaunay::init()
{
    setSelectRegionWidth(10);
    setSelectRegionHeight(10);
}
////////////////////////////////////////////////////////////////////////////////

EdgePtr Mesh:: addEdge( NodePtr &n0, NodePtr &n1, FacePtr &face)
{
    NodePtr vmin = std::min(n0,n1);

    for( auto oldedge : vmin->edges) {
        if( oldedge->hasNodes(n0,n1)) {
            oldedge->faces[1] = face;
            return oldedge;
        }
    }

    EdgePtr newedge(new Edge(n0,n1));
    newedge->faces[0] = face;
    vmin->edges.push_back(newedge);
    edges.push_back(newedge);
    return newedge;
}

////////////////////////////////////////////////////////////////////////////////

void Mesh:: makeConsistent( FacePtr &f)
{
    if( !f->active) return;

    int nCount = 0;

    NodePtr steiner[3];
    for( int i = 0; i < 3; i++) {
        if(f->edges[i]->steinerNode ) nCount++;
    }
    if( nCount == 0) return;

    NodePtr nodes[5];
    if( nCount == 1) {
        for( int i = 0; i < 3; i++) {
            if( f->edges[i]->steinerNode ) {
                nodes[0] = f->nodes[i];
                nodes[1] = f->nodes[(i+1)%3];
                nodes[2] = f->nodes[(i+2)%3];
                nodes[3] = f->edges[i]->steinerNode;
            }
        }
        remove(f);
        FacePtr f0 = Face::newObject(nodes[0], nodes[3], nodes[2]);
        FacePtr f1 = Face::newObject(nodes[1], nodes[2], nodes[3]);
        addFace(f0);
        addFace(f1);
        return;
    }

    if( nCount == 2) {
        for( int i = 0; i < 3; i++) {
            if( f->edges[i]->steinerNode == nullptr) {
                nodes[0] = f->nodes[i];
                nodes[1] = f->nodes[(i+1)%3];
                nodes[2] = f->nodes[(i+2)%3];
                nodes[3] = f->edges[(i+1)%3]->steinerNode;
                nodes[4] = f->edges[(i+2)%3]->steinerNode;
                break;
            }
        }
        remove(f);
        FacePtr f0 = Face::newObject( nodes[2], nodes[4], nodes[3]);
        FacePtr f1 = Face::newObject( nodes[0], nodes[1], nodes[3]);
        FacePtr f2 = Face::newObject( nodes[0], nodes[3], nodes[4]);
        addFace(f0);
        addFace(f1);
        addFace(f2);
        return;
    }
}

void IntrinsicDelaunay:: readMesh( const string &filename)
{
    std::vector<float>  nodes;
    std::vector<size_t> faces;

    ifstream ifile( filename.c_str(), ios::in);
    if( ifile.fail() ) {
        cout << "Warning: Input file not read " << endl;
        return;
    }
    string str;
    ifile >> str;
    if( str != "OFF") {
        cout << "Warning: Input file not in Off format" << endl;
        return;
    }

    size_t numNodes, numFaces, numEdges;
    ifile >> numNodes >> numFaces >> numEdges;

    double minval[3], maxval[3];

    minval[0] = std::numeric_limits<double>::max();
    minval[1] = std::numeric_limits<double>::max();
    minval[2] = std::numeric_limits<double>::max();

    maxval[0] =-std::numeric_limits<double>::max();
    maxval[1] =-std::numeric_limits<double>::max();
    maxval[2] =-std::numeric_limits<double>::max();

    double x, y, z;
    nodes.resize(3*numNodes);
    for( size_t i = 0; i < numNodes; i++) {
        ifile >> x >> y >> z;
        minval[0] = min( minval[0], x);
        minval[1] = min( minval[1], y);
        minval[2] = min( minval[2], z);

        maxval[0] = max( maxval[0], x);
        maxval[1] = max( maxval[1], y);
        maxval[2] = max( maxval[2], z);

        nodes[3*i+0] = x;
        nodes[3*i+1] = y;
        nodes[3*i+2] = z;
    }

    double xc = 0.5*(maxval[0] + minval[0]);
    double yc = 0.5*(maxval[1] + minval[1]);
    double zc = 0.5*(maxval[2] + minval[2]);

    double xlen = maxval[0] - minval[0];
    double ylen = maxval[1] - minval[1];
    double zlen = maxval[2] - minval[2];

    mesh.center[0] = xc;
    mesh.center[1] = yc;
    mesh.center[2] = zc;
    mesh.radius    = sqrt(xlen*xlen + ylen*ylen + zlen*zlen);

    mesh.nodes.resize(numNodes);
    for( size_t i = 0; i < numNodes; i++) {
        NodePtr v = Node::newObject();
        v->xyz[0] = nodes[3*i];
        v->xyz[1] = nodes[3*i+1];
        v->xyz[2] = nodes[3*i+2];
        v->id     = i;
        mesh.nodes[i] = v;
    }

    faces.resize(3*numFaces);
    size_t index = 0;
    int dummy, v0, v1, v2;

    for( size_t i = 0; i < numFaces; i++) {
        ifile >> dummy >> v0 >> v1 >> v2;
        assert( dummy == 3);
        faces[3*i+0] = v0;
        faces[3*i+1] = v1;
        faces[3*i+2] = v2;
        NodePtr n0   = mesh.nodes[v0];
        NodePtr n1   = mesh.nodes[v1];
        NodePtr n2   = mesh.nodes[v2];
        FacePtr newface = Face::newObject(n0,n1,n2);
        newface->id     = i;
        mesh.addFace(newface);
    }
}

////////////////////////////////////////////////////////////////////////////////

void Mesh::addFace( FacePtr &newface)
{
    assert( newface );
    auto n0 = newface->nodes[0]; assert( n0 );
    auto n1 = newface->nodes[1]; assert( n1 );
    auto n2 = newface->nodes[2]; assert( n2 );
    assert((n0 != n1) && (n1 != n2) && (n2 != n0));

    newface->edges[0] = addEdge(n0,n1,newface);
    newface->edges[1] = addEdge(n1,n2,newface);
    newface->edges[2] = addEdge(n2,n0,newface);

    n0->faces.push_back(newface);
    n1->faces.push_back(newface);
    n2->faces.push_back(newface);

    faces.push_back(newface);
}

////////////////////////////////////////////////////////////////////////////////
void Mesh::remove( FacePtr &oldface)
{
    for( int i = 0; i < 3; i++) {
        auto e = oldface->edges[i];
        if( e->faces[0] == oldface) {
            e->faces[0] = e->faces[1];
            e->faces[1] = nullptr;
        }
        if( e->faces[1] == oldface) {
            e->faces[1] = nullptr;
        }
    }

    for( int i = 0; i < 3; i++) {
        auto v  = oldface->nodes[i];
        if( !v->faces.empty() ) {
            auto it = std::remove(v->faces.begin(), v->faces.end(), oldface);
            v->faces.erase(it, v->faces.end() );
        }
    }

    oldface->active = 0;

}
////////////////////////////////////////////////////////////////////////////////
int Mesh::flip( EdgePtr &e)
{
    if( !e->active ) return 1;
    if( e->faces[0] == nullptr || e->faces[0] == nullptr) return 1;

    auto n0  = e->nodes[0];
    auto n1  = e->nodes[1];
    auto f0  = e->faces[0];
    auto f1  = e->faces[1];
    auto on0 = f0->getOpposite( n0, n1);
    auto on1 = f1->getOpposite( n0, n1);
    double d0 = length2(n0,n1);
    double d1 = length2(on0,on1);

    if( d0 < d1) return 1;

    double theta0 = f0->getAngleAt(on0);
    double theta1 = f1->getAngleAt(on1);

    if( theta0 + theta1 < 180) return 1;

    remove(f0);
    remove(f1);

    f0 = Face::newObject(n0, on0, on1);
    f1 = Face::newObject(n1, on1, on0);

    addFace(f0);
    addFace(f1);

    e->active = 0;
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
size_t Mesh::flip()
{
    size_t numEdgesFlipped = 0;
    while(1) {
        size_t numEdges = edges.size();
        size_t nCount = 0;
        for( size_t i = 0; i < numEdges; i++) {
            int err = flip( edges[i] );
            if( !err) nCount++;
        }
        numEdgesFlipped += nCount;
        if( nCount == 0) break;
    }
    return numEdgesFlipped;
}
////////////////////////////////////////////////////////////////////////////////

void Mesh::refine(FacePtr &f, int type)
{
    assert(f);

    if(f->active) {
        if(f->getArea() > 1.0E-10) {
            remove(f);

            if( type == 14) {
                NodePtr steiner[3];
                for( int i = 0; i < 3; i++) {
                    auto edge = f->edges[i];
                    if( edge->steinerNode == nullptr) {
                        auto newnode = Node::newObject();
                        newnode->xyz = edge->getCenter();
                        newnode->id  = nodes.size();
                        nodes.push_back(newnode);
                        edge->steinerNode = newnode;
                    }
                    steiner[i] = edge->steinerNode;
                    assert( steiner[i] );
                }
                FacePtr f0 = Face::newObject( f->nodes[0], steiner[0], steiner[2]);
                FacePtr f1 = Face::newObject( f->nodes[1], steiner[1], steiner[0]);
                FacePtr f2 = Face::newObject( f->nodes[2], steiner[2], steiner[1]);
                FacePtr f3 = Face::newObject( steiner[0],  steiner[1], steiner[2]);
                addFace(f0);
                addFace(f1);
                addFace(f2);
                addFace(f3);
            }

            if( type == 13) {
                NodePtr newnode = Node::newObject();
                newnode->xyz    = f->getCentroid();
                newnode->id     = nodes.size();
                nodes.push_back(newnode);

                FacePtr f0 = Face::newObject( f->nodes[0], f->nodes[1], newnode);
                FacePtr f1 = Face::newObject( f->nodes[1], f->nodes[2], newnode);
                FacePtr f2 = Face::newObject( f->nodes[2], f->nodes[0], newnode);
                addFace(f0);
                addFace(f1);
                addFace(f2);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
void Mesh::refine(int type)
{
    size_t numFaces = faces.size();
    for( size_t i = 0; i < numFaces; i++)
        refine( faces[i], type );
}
////////////////////////////////////////////////////////////////////////////////

void Mesh::saveAs( const std::string &filename)
{
    ofstream ofile(filename.c_str(), ios::out);

    size_t  numnodes = 0;
    for( auto v: nodes) {
        if( v->active ) v->id = numnodes++;
    }

    set<NodePtr> vSet;
    size_t  numfaces = 0;
    for( auto f: faces) {
        if( f->active ) {
            vSet.insert(f->nodes[0]);
            vSet.insert(f->nodes[1]);
            vSet.insert(f->nodes[2]);
            numfaces++;
        }
    }

    size_t index = 0;
    for( auto v: vSet) {
        v->id = index++;
        ofile << "v " << v->xyz[0] << " " << v->xyz[1] << " " << v->xyz[2] << endl;
    }

    for( auto f: faces) {
        if( f->active )
            ofile << "f " << f->nodes[0]->id +1 << " "
                  << f->nodes[1]->id +1 << " "
                  << f->nodes[2]->id +1 << endl;
    }

    cout << "Info: Cutout written to file: " << filename << endl;
}

////////////////////////////////////////////////////////////////////////////////

void IntrinsicDelaunay::keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_0) {
        pickEntity = 0;
        this->setSelectedName(-1);
    }

    if( e->key() == Qt::Key_R) {
        mesh.refine(13);
        update();
        return;
    }

    if( e->key() == Qt::Key_X) {
        mesh.flip();
        update();
        return;
    }

    if( e->key() == Qt::Key_F) {
        mesh.saveAs("cutout.obj");
        update();
        return;
    }

    if( e->key() == Qt::Key_N) {
        displayIDs = !displayIDs;
    }

    if( e->key() == Qt::Key_S) {
        displaySurface = !displaySurface;
        update();
        return;
    }

    if( e->key() == Qt::Key_W) {
        displayWires = !displayWires;
    }

    if( e->key() == Qt::Key_L) {
        useLights = !useLights;
        update();
        return;
    }


    if( e->key() == Qt::Key_Home) {
        qglviewer::Vec pos;
        pos[0]  = mesh.center[0];
        pos[1]  = mesh.center[1];
        pos[2]  = mesh.center[2];
        camera()->setSceneCenter(pos);
        camera()->setSceneRadius(mesh.radius);
        camera()->centerScene();
        camera()->showEntireScene();
        update();
        return;
    }

    QGLViewer::keyPressEvent(e);

    update();
}

////////////////////////////////////////////////////////////////////////////////

void IntrinsicDelaunay:: mousePressEvent( QMouseEvent *e)
{
    QGLViewer::mousePressEvent(e);
}

////////////////////////////////////////////////////////////////////////////////

void IntrinsicDelaunay:: selectNode( size_t id)
{
}
////////////////////////////////////////////////////////////////////////////////

void IntrinsicDelaunay::mouseReleaseEvent( QMouseEvent *e)
{
    int id = this->selectedName();

    if( id >= 0) selectNode(id);

    QGLViewer::mouseReleaseEvent(e);

    update();
}

////////////////////////////////////////////////////////////////////////////////
void IntrinsicDelaunay::drawNodes()
{
    glDisable( GL_LIGHTING);
    glPointSize(2);
    glColor3f( 0.0, 0.0, 1.0);

    /*
    glBegin(GL_POINTS);
    for( auto v: mesh.nodes) {
        if(v->active) glVertex3fv( &v->xyz[0]);
    }
    glEnd();

    if( !contours.empty() ) {
        glPointSize(5);
        glColor3f( 1.0, 0.0, 0.0);
        glBegin(GL_POINTS);
        for ( size_t i = 0; i < contours.size(); i++) {
            for ( auto v: contours[i].boundnodes) {
                glVertex3fv( &v->xyz[0] );
            }
        }
        glEnd();
    }

    if( !newContour.boundnodes.empty() ) {
        glPointSize(5);
        glColor3f( 1.0, 0.0, 0.0);
        glBegin(GL_POINTS);
        for ( auto v: newContour.boundnodes) {
            glVertex3fv( &v->xyz[0] );
        }
        glEnd();
    }

    glPointSize(5);
    glColor3f( 1.0, 1.0, 0.0);
    glBegin(GL_POINTS);
    for ( size_t j = 0; j < contours.size(); j++) {
        for( auto v: contours[j].landmarks)
            glVertex3fv( &v->xyz[0] );
    }
    */
}

////////////////////////////////////////////////////////////////////////////////

void IntrinsicDelaunay::drawFaces(int filled)
{
    if( useLights ) glEnable(GL_LIGHTING);

    if( filled ) {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glColor3f( 0.8, 0.9, 0.9);
        for ( auto f : mesh.faces) {
            if(f->active) {
                glBegin(GL_TRIANGLES);
                glVertex3fv( &f->nodes[0]->xyz[0] );
                glVertex3fv( &f->nodes[1]->xyz[0] );
                glVertex3fv( &f->nodes[2]->xyz[0] );
                glEnd();
            }
        }
    } else {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glColor3f( 0.2, 0.2, 0.2);
        for ( auto f : mesh.faces) {
            if(f->active) {
                glBegin(GL_LINE_LOOP);
                glVertex3fv( &f->nodes[0]->xyz[0] );
                glVertex3fv( &f->nodes[1]->xyz[0] );
                glVertex3fv( &f->nodes[2]->xyz[0] );
                glEnd();
            }
        }
    }



    if( pickEntity == 2) {
        int id = this->selectedName();
        if( id >= 0) {
            auto f = mesh.faces[id];
            if(f->active) {
                glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                glColor3f( 0.0, 1.0, 0.0);
                glBegin(GL_TRIANGLES);
                glVertex3fv( &f->nodes[0]->xyz[0] );
                glVertex3fv( &f->nodes[1]->xyz[0] );
                glVertex3fv( &f->nodes[2]->xyz[0] );
                glEnd();
                glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
void IntrinsicDelaunay:: drawIDs()
{
}
////////////////////////////////////////////////////////////////////////////////

void IntrinsicDelaunay::drawWithNames()
{
    glDisable( GL_LIGHTING);
    glPointSize(2);
    glColor3f( 0.0, 0.0, 1.0);

    if( pickEntity == 0) {
        for ( auto v: mesh.nodes) {
            if( v->active) {
                glPushName(v->id);
                glBegin(GL_POINTS);
                glVertex3fv( &v->xyz[0] );
                glEnd();
                glPopName();
            }
        }
    }

    if( pickEntity == 2) {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glColor3f( 0.8, 0.9, 0.9);
        for ( auto f: mesh.faces) {
            if( f->active) {
                glPushName(f->id);
                glBegin(GL_TRIANGLES);
                glVertex3fv( &f->nodes[0]->xyz[0] );
                glVertex3fv( &f->nodes[1]->xyz[0] );
                glVertex3fv( &f->nodes[2]->xyz[0] );
                glEnd();
                glPopName();
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void IntrinsicDelaunay::draw()
{
    glPolygonOffset(1.0,1.0);
    glEnable(GL_POLYGON_OFFSET_FILL);

    drawNodes();
    if( displayWires)   drawFaces(0);
    if( displaySurface) drawFaces(1);
    if( displayIDs)     drawIDs();
}

////////////////////////////////////////////////////////////////////////////////
