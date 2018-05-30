#ifndef SRC_FAMILY_H
#define SRC_FAMILY_H

#include "utils.h"

class Family {
private:
    FamilyType family_type;
    LinkType link_type;
public:
    Family(FamilyType family_type);

    Family(FamilyType family_type, LinkType link_type);

    vec (*link_fun)(const vec &mu);

    vec (*link_inv)(const vec &eta);

    vec (*derivative)(const vec &eta);

    vec (*variance)(const vec &mu);

    vec (*likelyhood)(const vec &y, const vec &mu);
};


#endif //SRC_FAMILY_H
